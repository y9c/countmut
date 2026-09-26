/* countmut_core.c -- C computation core for the unified countmut tool.
 *
 * Implements a position-by-position pileup walk (htslib bam_mplp_auto) and
 * counts, per genomic site and per biological strand, the observed bases binned
 * by quality/conversion category (mutation mode) or by residue (base/allele
 * mode), with mate-overlap deduplication by query name and parallel processing
 * across genomic bins.
 *
 * Semantics distilled from:
 *   - minipileup (pileup walk, read + base filters, allele counting)
 *   - perbase / pbr (mate-aware overlap dedup, PileupPosition a/c/g/t/n/ins/del)
 *   - countmut (biological strand, trim orientation, bisulfite NS/Zf/Yf tiers)
 *
 * Author: Ye Chang
 * Date: 2026-08-30
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <assert.h>
#include <pthread.h>
#include <unistd.h>
#include <zlib.h>

#include <htslib/sam.h>
#ifndef bam_pe32
static inline uint32_t bam_pe32(const void *p) { uint32_t v; memcpy(&v, p, 4); return v; }
#endif
#ifndef bam_cigar_op_p
#define bam_cigar_op_p(p) bam_cigar_op(bam_pe32(p))
#endif
#ifndef bam_cigar_oplen_p
#define bam_cigar_oplen_p(p) bam_cigar_oplen(bam_pe32(p))
#endif
#include <htslib/faidx.h>
#include <htslib/thread_pool.h>
#include "ksort.h"
#include "khash.h"
#include "countmut_core.h"
#include "countmut_expr.h"

/* local flag constants not always exposed by the bundled sam.h */
#ifndef BAM_FPAIRED
#define BAM_FPAIRED 1
#endif
#ifndef BAM_FQCFAIL
#define BAM_FQCFAIL 512
#endif

/* nt16 nucleotide codes -> base index (0=A,1=C,2=G,3=T,4=N) */
static int nt16_index(uint8_t b) {
    switch (b) {
    case 1: return 0;  /* A */
    case 2: return 1;  /* C */
    case 4: return 2;  /* G */
    case 8: return 3;  /* T */
    case 15: return 4; /* N */
    default: return 4;
    }
}

/* biological strand: 0 = '+', 1 = '-' */
static int bio_strand(const bam1_t *b) {
    int rev = bam_is_rev(b);
    if (b->core.flag & BAM_FPAIRED) {
        if (b->core.flag & BAM_FREAD1) return rev ? 1 : 0;
        return rev ? 0 : 1;
    }
    return rev ? 1 : 0;
}

/* Read (R1/R2) query-end trimming, in addition to the fragment 5'/3' trim.
 *   R1: `r1_end` bases off its 3' query end  (qpos >= len - r1_end)
 *   R2: `r2_start` bases off its 5' query start (qpos < r2_start)
 * Returns 1 when the base at `qpos` should be skipped. */
/* per-site accumulator: site_t lives in countmut_core.h (cnt is
 * [strand][CM_CAT_MAX][base]) so the expression layer can index it too. */
static void site_zero(site_t *s) { memset(s, 0, sizeof(*s)); }

static int better(int mapq, int r1, int q, int omapq, int or1, int oq) {
    if (mapq != omapq) return mapq > omapq;
    if (r1 != or1) return r1 > or1;
    return q > oq;
}

static int base_to_index(char c) {
    switch (toupper((unsigned char)c)) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return 4;
    }
}

/* Reverse-complement a single uppercase motif base (A<->T, C<->G, else N). */
static char rc_nt(char c) {
    switch (c) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    default: return 'N';
    }
}

/* read-level filters (samtools reqflags/exclflags; mapq/quality are -e filters) */
static int read_fails(const cm_config *cfg, const bam1_t *b) {
    if (cfg->req_flags && (b->core.flag & (uint32_t)cfg->req_flags) != (uint32_t)cfg->req_flags) return 1;
    if (cfg->excl_flags && (b->core.flag & (uint32_t)cfg->excl_flags)) return 1;
    return 0;
}

typedef struct {
    samFile *fp;
    sam_hdr_t *hdr;
    hts_itr_t *itr;
    int beg, end;
} aux_t;

static int read_bam(void *data, bam1_t *b) {
    aux_t *aux = (aux_t *)data;
    int ret = aux->itr ? sam_itr_next(aux->fp, aux->itr, b) : sam_read1(aux->fp, aux->hdr, b);
    return ret;
}

KHASH_INIT(qn, char *, int, 1, kh_str_hash_func, kh_str_hash_equal)

/* -e read-constant memo for the pileup engine, keyed by the pileup slot
 * pointer (p->b): htslib keeps the same bam1_t for one read across all the
 * positions it covers, so an int-keyed slot hash lets us evaluate a
 * read-constant expression once per read and reuse it for every appearance.
 * A pos/qlen/qname verify guards against a recycled buffer now holding a
 * different read.  Each mplp slot holds at most one read at a time. */
typedef struct { int64_t pos; int qlen; char *qn; int slot; } expr_cc_t;
static inline khint_t pex_hash(uintptr_t p) {
    return (khint_t)(p >> 3) ^ ((khint_t)(p >> 13) & 0x0ff);
}
static inline int pex_equal(uintptr_t a, uintptr_t b) { return a == b; }
KHASH_INIT(pex, uintptr_t, expr_cc_t, 1, pex_hash, pex_equal)

/* BED / position-list region support (from bedidx.c, mirrors minipileup) */
void *bed_read(const char *fn);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void bed_destroy(void *_h);

/* per-worker reusable state */
typedef struct {
    khash_t(qn) *kh;
    int *sel, *mapq_a, *r1_a, *q_a, *g_a, sel_cap;
    char *motif_buf;      /* reference-forward motif window (per-site) */
    char *motif_rc_buf;   /* reverse-complemented copy for the minus-strand row */
    char *chr_seq; int chr_len, last_tid;
    samFile *fp; hts_idx_t *idx; faidx_t *fai;
    htsThreadPool *tpool;   /* shared CRAM decode pool (NULL = no pool) */
    void *inc_bed, *exc_bed;
    cm_expr *expr;   /* Lua -e / -p filters (NULL when none) */
    khash_t(pex) *pexc; /* -e read-constant memo keyed by pileup slot (pileup) */
} worker_t;

static void worker_init(worker_t *w, const char *bam, const char *fa, int pad,
                        const char *bedfile, const char *exclude,
                        const char *read_expr, const char *pile_expr,
                        const char *output_expr, htsThreadPool *tpool) {
    w->kh = kh_init(qn);
    w->sel = w->mapq_a = w->r1_a = w->q_a = w->g_a = NULL;
    w->sel_cap = 0;
    w->motif_buf = (char *)malloc(2 * pad + 2);
    w->motif_rc_buf = (char *)malloc(2 * pad + 2);
    w->chr_seq = NULL; w->chr_len = 0; w->last_tid = -1;
    w->fp = sam_open(bam, "r");
    hts_set_fai_filename(w->fp, fa);
    if (tpool) hts_set_thread_pool(w->fp, tpool);   /* shared CRAM decode pool */
    hts_set_threads(w->fp, 2);   /* CRAM decode needs the reference */
    w->idx = sam_index_load(w->fp, bam);
    w->fai = fai_load(fa);
    w->inc_bed = bedfile ? bed_read(bedfile) : NULL;
    w->exc_bed = exclude ? bed_read(exclude) : NULL;
    w->expr = cm_expr_new(read_expr, pile_expr, output_expr);
    w->pexc = kh_init(pex);
}

static void worker_free(worker_t *w) {
    if (w->sel) free(w->sel);
    if (w->mapq_a) free(w->mapq_a);
    if (w->r1_a) free(w->r1_a);
    if (w->q_a) free(w->q_a);
    if (w->g_a) free(w->g_a);
    if (w->motif_buf) free(w->motif_buf);
    if (w->motif_rc_buf) free(w->motif_rc_buf);
    if (w->chr_seq) free(w->chr_seq);
    kh_destroy(qn, w->kh);
    if (w->exc_bed) bed_destroy(w->exc_bed);
    if (w->inc_bed) bed_destroy(w->inc_bed);
    if (w->fai) fai_destroy(w->fai);
    if (w->idx) hts_idx_destroy(w->idx);
    if (w->fp) sam_close(w->fp);
    cm_expr_free(w->expr);
    if (w->pexc) {
        for (khint_t k = kh_begin(w->pexc); k != kh_end(w->pexc); ++k)
            if (kh_exist(w->pexc, k)) free(kh_val(w->pexc, k).qn);
        kh_destroy(pex, w->pexc);
    }
}

/* -e read filter, memoized by (ref_start, qname) WHEN the expression is
 * read-constant (no per-base qpos/bq/base/ref/dist).  A read-constant
 * expression yields the same result at every base of a read, so caching makes
 * the pileup engine evaluate it once per read instead of once per position
 * (the read-walk engine already does once per read).  Per-base expressions are
 * evaluated at every aligned base, uncached. */
/* Route one aligned base: returns the category slot to count it into
 * (-1 = drop, 0..CM_CAT_MAX-1 = slot).  The router subsumes both the old -e
 * keep/drop decision (nil/false -> -1) and the per-category tier assignment
 * (a number -> that slot; true -> default slot 0). */
static int expr_pass(worker_t *w, const bam1_t *b, const char *rname,
                     const char *mrname, int s, int qpos, char ref_ch) {
    cm_expr *x = w->expr;
    if (x == NULL || !cm_expr_has_read(x)) return 0;   /* no -e -> default slot 0 */
    if (!cm_expr_read_constant(x))
        return cm_expr_route(x, b, rname, mrname, qpos, s ? -1 : 1, ref_ch);
    /* read-constant: memoize by the pileup slot pointer (stable per read
     * across its span) with a pos/qlen/qname verify against recycling. */
    uintptr_t key = (uintptr_t)(const void *)b;
    khint_t k = kh_get(pex, w->pexc, key);
    if (k != kh_end(w->pexc)) {
        expr_cc_t *cc = &kh_val(w->pexc, k);
        if (cc->pos == b->core.pos && cc->qlen == (int)b->core.l_qseq
            && cc->qn && strcmp(cc->qn, bam_get_qname(b)) == 0)
            return cc->slot;
        free(cc->qn); cc->qn = NULL;
    }
    int slot = cm_expr_route(x, b, rname, mrname, 0, s ? -1 : 1, 'N');
    if (k == kh_end(w->pexc)) {
        int ret; k = kh_put(pex, w->pexc, key, &ret);
        memset(&kh_val(w->pexc, k), 0, sizeof(expr_cc_t)); /* fresh slots are uninitialized */
    }
    expr_cc_t *cc = &kh_val(w->pexc, k);
    if (cc->qn == NULL) cc->qn = strdup(bam_get_qname(b));
    cc->pos = b->core.pos;
    cc->qlen = (int)b->core.l_qseq;
    cc->slot = slot;
    return slot;
}

typedef struct { int tid, beg, end; } region_t;

/* record of which worker owns which region's rows in that worker's temp file,
 * so the final output can be re-assembled in global region order regardless of
 * the (dynamic) claim order */
typedef struct { int worker; long off; long len; } span_t;

typedef struct {
    const cm_config *cfg;
    bam_hdr_t *hdr;
    const char *bam;
    region_t *regions;
    int nregions;
    volatile int done;
    volatile int next;   /* next region to claim (dynamic work queue) */
    span_t *spans;
    FILE **files;
    worker_t *workers;
    int nthreads;
} work_t;

/* Write the header line for the configured output format. */
static void write_header(FILE *fp, const cm_config *cfg) {
    if (cfg->output_expr && cfg->output_expr[0]) {   /* custom -o template */
        if (cfg->fmt_header && cfg->fmt_header[0]) {
            for (const char *p = cfg->fmt_header; *p; ++p) {   /* expand \t \n \\ */
                if (*p == '\\' && p[1] == 't') { fputc('\t', fp); ++p; }
                else if (*p == '\\' && p[1] == 'n') { fputc('\n', fp); ++p; }
                else if (*p == '\\' && p[1] == '\\') { fputc('\\', fp); ++p; }
                else fputc(*p, fp);
            }
            fputc('\n', fp);
        }
        return;
    }
    if (cfg->out == CM_OUT_COMPOSITION) {
        if (cfg->strandless) fputs("chrom\tpos\tref\tdepth\ta\tc\tg\tt\tn", fp);
        else fputs("chrom\tpos\tstrand\tref\tdepth\ta\tc\tg\tt\tn", fp);
        if (cfg->count_indels) fputs("\tins\tdel\tref_skip\tfail", fp);
        fputc('\n', fp);
    } else {
        if (cfg->vcf) {
            fputs("##fileformat=VCFv4.2\n", fp);
            fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n", fp);
        } else {
            fputs("chrom\tpos\tref\tdepth\tref_count\talt\talt_count\n", fp);
        }
    }
}

/* Emit one site -- shared by the pileup engine and the read-walk engine so the
 * two walks produce identical rows. */
static void emit_site(worker_t *w, const cm_config *cfg, bam_hdr_t *hdr, FILE *fp,
                      int tid, int64_t pos, char ref_ch, const site_t *site,
                      int emit_plus, int emit_minus) {
    /* custom output template (-o): evaluate it per emitted strand and write
     * exactly what it returns; bypasses the built-in formats. */
    if (w != NULL && w->expr && cm_expr_has_output(w->expr)) {
        /* Reference-forward motif window for {motif} (built whenever the
         * sequence is available, not only in the legacy conversion view);
         * the minus-strand row gets its reverse complement (parity with the
         * bisulfite reformat bridge). */
        const char *motif = NULL;
        if (w->chr_len > 0) {
            int mlen = cfg->pad * 2 + 1;
            for (int k2 = (int)pos - cfg->pad; k2 < (int)pos + cfg->pad + 1; ++k2)
                w->motif_buf[k2 - ((int)pos - cfg->pad)] =
                    (k2 < 0 || k2 >= w->chr_len) ? 'N' : (char)toupper((unsigned char)w->chr_seq[k2]);
            w->motif_buf[mlen] = 0;
            motif = w->motif_buf;
        }
        /* No ref/mut targets in the template path (unified counter): the
         * derived u/m/o fields stay off for template rows. */
        int refi = -1, muti = -1;
        if (cfg->strandless) {
            /* Strandless: combine both strands into one row (parity with the
             * built-in composition format).  The motif is reference-forward
             * (no reverse complement needed). */
            int cnt[CM_CAT_MAX][5] = {0};
            int ins = 0, del = 0, rs = 0, fl = 0;
            if (emit_plus) {
                for (int c = 0; c < CM_CAT_MAX; ++c)
                    for (int b = 0; b < 5; ++b) cnt[c][b] += site->cnt[0][c][b];
                ins += site->ins[0]; del += site->del[0]; rs += site->refskip[0]; fl += site->fail[0];
            }
            if (emit_minus) {
                for (int c = 0; c < CM_CAT_MAX; ++c)
                    for (int b = 0; b < 5; ++b) cnt[c][b] += site->cnt[1][c][b];
                ins += site->ins[1]; del += site->del[1]; rs += site->refskip[1]; fl += site->fail[1];
            }
            int sdepth = ins + del + rs + fl;
            for (int c = 0; c < CM_CAT_MAX; ++c)
                for (int b = 0; b < 5; ++b) sdepth += cnt[c][b];
            if (sdepth == 0) return;
            cm_expr_output(w->expr, hdr->target_name[tid], pos, ref_ch, motif,
                           cnt, ins, del, rs, fl, refi, muti, -1, fp);
        } else {
            /* Per-strand rows: counts are this strand only, per category slot. */
            for (int s = 0; s < 2; ++s) {
                if (s == 0 && !emit_plus) continue;
                if (s == 1 && !emit_minus) continue;
                /* target_base: emit only the strand whose reference base (on that
                 * strand) equals the target.  For m6A (target A) this keeps the
                 * '+' strand A-site (genomic ref A) and the '-' strand A-site
                 * (genomic ref T), dropping the spurious complement-strand rows. */
                if (cfg->target_base >= 0) {
                    char sref = (s == 0) ? ref_ch : rc_nt(ref_ch);
                    if (base_to_index((char)sref) != cfg->target_base) continue;
                }
                int sdepth = site->ins[s] + site->del[s] + site->refskip[s] + site->fail[s];
                for (int c = 0; c < CM_CAT_MAX; ++c)
                    for (int b = 0; b < 5; ++b) sdepth += site->cnt[s][c][b];
                if (sdepth == 0) continue;   /* match the built-in formats: skip empty strands */
                const char *mot_s = motif;
                if (motif && s == 1) {       /* minus strand: reverse complement */
                    int mlen = cfg->pad * 2 + 1;
                    for (int k2 = 0; k2 < mlen; ++k2)
                        w->motif_rc_buf[k2] = rc_nt(motif[mlen - 1 - k2]);
                    w->motif_rc_buf[mlen] = 0;
                    mot_s = w->motif_rc_buf;
                }
                cm_expr_output(w->expr, hdr->target_name[tid], pos, ref_ch, mot_s,
                               site->cnt[s], site->ins[s], site->del[s],
                               site->refskip[s], site->fail[s], refi, muti, s, fp);
            }
        }
        return;
    }
    if (cfg->out == CM_OUT_COMPOSITION) {
        if (!cfg->strandless) {
            for (int s = 0; s < 2; ++s) {
                if (s == 0 && !emit_plus) continue;
                if (s == 1 && !emit_minus) continue;
                int cnt[5] = {0};
                for (int c = 0; c < CM_CAT_MAX; ++c)
                    for (int b = 0; b < 5; ++b) cnt[b] += site->cnt[s][c][b];
                int dep = cnt[0]+cnt[1]+cnt[2]+cnt[3]+cnt[4];
                int t_ins = site->ins[s], t_del = site->del[s], t_rs = site->refskip[s], t_fl = site->fail[s];
                if (dep + t_rs + t_del + t_ins + t_fl == 0) continue;
                fprintf(fp, "%s\t%d\t%c\t%c\t%d\t%d\t%d\t%d\t%d\t%d",
                        hdr->target_name[tid], (int)pos + 1, s ? '-' : '+', ref_ch, dep,
                        cnt[0], cnt[1], cnt[2], cnt[3], cnt[4]);
                if (cfg->count_indels) fprintf(fp, "\t%d\t%d\t%d\t%d", t_ins, t_del, t_rs, t_fl);
                fputc('\n', fp);
            }
        } else {
            int cnt[5] = {0}; int t_ins = 0, t_del = 0, t_rs = 0, t_fl = 0;
            if (emit_plus) {
                for (int c = 0; c < CM_CAT_MAX; ++c)
                    for (int b = 0; b < 5; ++b) cnt[b] += site->cnt[0][c][b];
                t_ins += site->ins[0]; t_del += site->del[0]; t_rs += site->refskip[0]; t_fl += site->fail[0];
            }
            if (emit_minus) {
                for (int c = 0; c < CM_CAT_MAX; ++c)
                    for (int b = 0; b < 5; ++b) cnt[b] += site->cnt[1][c][b];
                t_ins += site->ins[1]; t_del += site->del[1]; t_rs += site->refskip[1]; t_fl += site->fail[1];
            }
            int dep = cnt[0]+cnt[1]+cnt[2]+cnt[3]+cnt[4];
            if (dep + t_rs + t_del + t_ins + t_fl == 0) return;
            fprintf(fp, "%s\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d",
                    hdr->target_name[tid], (int)pos + 1, ref_ch, dep, cnt[0], cnt[1], cnt[2], cnt[3], cnt[4]);
            if (cfg->count_indels) fprintf(fp, "\t%d\t%d\t%d\t%d", t_ins, t_del, t_rs, t_fl);
            fputc('\n', fp);
        }
    } else { /* allele */
        int cnt[5] = {0};
        if (emit_plus)
            for (int c = 0; c < CM_CAT_MAX; ++c)
                for (int b = 0; b < 5; ++b) cnt[b] += site->cnt[0][c][b];
        if (emit_minus)
            for (int c = 0; c < CM_CAT_MAX; ++c)
                for (int b = 0; b < 5; ++b) cnt[b] += site->cnt[1][c][b];
        int dep = cnt[0]+cnt[1]+cnt[2]+cnt[3]+cnt[4];
        if (dep <= 0) return;
        int refi = base_to_index(ref_ch), ref_n = cnt[refi], best = -1, bn = 0;
        for (int i = 0; i < 4; ++i) if (i != refi && cnt[i] > bn) { bn = cnt[i]; best = i; }
        if (cfg->vcf) {
            if (best < 0) return;
            const char *alts = "ACGT";
            fprintf(fp, "%s\t%d\t.\t%c\t%c\t.\tPASS\t.\tGT:AD\t0/1:%d,%d\n",
                    hdr->target_name[tid], (int)pos + 1, ref_ch, alts[best], ref_n, bn);
        } else {
            fprintf(fp, "%s\t%d\t%c\t%d\t%d\t%c\t%d\n",
                    hdr->target_name[tid], (int)pos + 1, ref_ch, dep, ref_n,
                    best < 0 ? '.' : "ACGT"[best], best < 0 ? 0 : bn);
        }
    }
}

/* Evaluate the -p (pile/site) filter for a fully-built site_t.  Computes the
 * A/C/G/T/N totals (both strands, all quality tiers) plus indels, and the
 * reference window for mutation mode.  Returns 1 = keep, 0 = omit. */
/* Evaluate the -p site filter PER STRAND.  Returns a bitmask: bit 0 = the
 * '+' strand passes, bit 1 = the '-' strand passes.  A site is kept if at
 * least one strand passes; emit_site then emits only the passing strand(s),
 * so a strand-aware filter like `base == 'A'` keeps only the A-site strand
 * and drops the spurious complement-strand rows.  When -p is not set both
 * bits are set (no strand filtering). */
static int expr_pile_apply_strands(cm_expr *x, const cm_config *cfg, worker_t *w,
                                   const char *chrom, const site_t *site,
                                   int64_t pos, char ref_ch) {
    if (x == NULL || !cm_expr_has_pile(x)) return 3;   /* both strands */
    const char *motif = NULL;
    if (w->chr_len > 0) {   /* motif window for -p/-o expressions (ref-forward) */
        int mlen = cfg->pad * 2 + 1;
        for (int k2 = (int)pos - cfg->pad; k2 < (int)pos + cfg->pad + 1; ++k2) {
            w->motif_buf[k2 - ((int)pos - cfg->pad)] =
                (k2 < 0 || k2 >= w->chr_len) ? 'N' : (char)toupper((unsigned char)w->chr_seq[k2]);
        }
        w->motif_buf[mlen] = 0;
        motif = w->motif_buf;
    }
    int refi = -1, muti = -1;   /* no ref/mut targets in the unified counter */
    int mask = 0;
    for (int s = 0; s < 2; ++s) {
        int cnt[5] = {0};
        int ins = site->ins[s], del = site->del[s], rs = site->refskip[s], fl = site->fail[s];
        for (int c = 0; c < CM_CAT_MAX; ++c)
            for (int b = 0; b < 5; ++b) cnt[b] += site->cnt[s][c][b];
        if (cm_expr_pile_strand(x, chrom, pos, ref_ch, motif, cnt, ins, del, rs, fl,
                                refi, muti, s))
            mask |= (1 << s);
    }
    return mask;
}

/* Fetch + uppercase the chromosome sequence once per tid (instead of calling
 * toupper() on every per-base / per-position access).  Output is unchanged. */
static void load_chr_seq(worker_t *w, bam_hdr_t *hdr, int tid) {
    if (w->chr_seq) free(w->chr_seq);
    w->chr_seq = fai_fetch(w->fai, hdr->target_name[tid], &w->chr_len);
    w->last_tid = tid;
    if (w->chr_seq)
        for (int i = 0; i < w->chr_len; ++i)
            w->chr_seq[i] = (char)toupper((unsigned char)w->chr_seq[i]);
}

/* Count one interval [beg,end) of `tid` and write rows to fp. */
static void count_interval(worker_t *w, const cm_config *cfg, bam_hdr_t *hdr, FILE *fp, int tid, int beg, int end) {
    aux_t aux;
    aux.fp = w->fp; aux.hdr = hdr;
    aux.beg = beg; aux.end = end;
    aux.itr = w->idx ? sam_itr_queryi(w->idx, tid, beg, end) : NULL;

    int *n_plp = (int *)calloc(1, sizeof(int));
    const bam_pileup1_t **plp = (const bam_pileup1_t **)calloc(1, sizeof(void *));
    void *data_ptrs[1] = {&aux};
    bam_mplp_t mplp = bam_mplp_init(1, read_bam, data_ptrs);
    /* htslib's default maxcnt is 8000; 0 = unlimited so we raise it. */
    bam_mplp_set_maxcnt(mplp, cfg->max_depth > 0 ? cfg->max_depth : 0x7fffffff);

    int pos;
    site_t site;
    while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) {
        if (pos < beg || pos >= end) continue;
        if (w->last_tid != tid) load_chr_seq(w, hdr, tid);
        if (w->chr_len == 0 || pos >= w->chr_len) continue;
        /* BED region restriction (pbr -b include / -x exclude) */
        if (w->inc_bed && !bed_overlap(w->inc_bed, hdr->target_name[tid], pos, pos + 1)) continue;
        if (w->exc_bed && bed_overlap(w->exc_bed, hdr->target_name[tid], pos, pos + 1)) continue;
        int n = n_plp[0];
        if (n == 0) continue;
        const char ref_ch = w->chr_seq[pos];   /* pre-uppercased */

        if (n > w->sel_cap) {
            w->sel = (int *)realloc(w->sel, n * sizeof(int));
            w->mapq_a = (int *)realloc(w->mapq_a, n * sizeof(int));
            w->r1_a = (int *)realloc(w->r1_a, n * sizeof(int));
            w->q_a = (int *)realloc(w->q_a, n * sizeof(int));
            w->g_a = (int *)realloc(w->g_a, n * sizeof(int));
            w->sel_cap = n;
        }
        site_zero(&site);
        memset(w->sel, 0, n * sizeof(int));
        kh_clear(qn, w->kh);
        for (int i = 0; i < n; ++i) {
            const bam_pileup1_t *p = &plp[0][i];
            const bam1_t *b = p->b;
            int s = bio_strand(b);
            if (cfg->strand_process == CM_STRAND_FORWARD && s != 0) continue;
            if (cfg->strand_process == CM_STRAND_REVERSE && s != 1) continue;
            if (read_fails(cfg, b)) { site.fail[s]++; continue; }
            if (p->is_refskip) { site.refskip[s]++; continue; }
            if (p->is_del) { site.del[s]++; continue; }
            if (p->qpos < 0 || p->qpos >= b->core.l_qseq) continue;
            int qpos = p->qpos;
            /* -e read router: returns the category slot (-1 = drop).  Evaluated
             * once per read when read-constant via the exprc memo, else per
             * aligned base (the same spot as the Python engine). */
            int g = expr_pass(w, b, hdr->target_name[tid],
                              (b->core.mtid >= 0 && b->core.mtid < hdr->n_targets)
                                  ? hdr->target_name[b->core.mtid] : "",
                              s, qpos, ref_ch);
            if (g < 0) continue;
            w->g_a[i] = g;
            int mapq = (int)b->core.qual;
            int r1 = (b->core.flag & BAM_FREAD1) ? 1 : 0;
            int qual = (int)bam_get_qual(b)[qpos];
            const char *qname = bam_get_qname(b);
            khint_t k = kh_get(qn, w->kh, qname);
            if (k == kh_end(w->kh)) {
                int ret; k = kh_put(qn, w->kh, (char *)qname, &ret);
                kh_val(w->kh, k) = i;
                w->sel[i] = 1; w->mapq_a[i] = mapq; w->r1_a[i] = r1; w->q_a[i] = qual;
            } else {
                int j = kh_val(w->kh, k);
                if (better(mapq, r1, qual, w->mapq_a[j], w->r1_a[j], w->q_a[j])) {
                    w->sel[j] = 0;
                    kh_val(w->kh, k) = i;
                    w->sel[i] = 1; w->mapq_a[i] = mapq; w->r1_a[i] = r1; w->q_a[i] = qual;
                } else {
                    w->sel[i] = 0; w->mapq_a[i] = mapq; w->r1_a[i] = r1; w->q_a[i] = qual;
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            if (!w->sel[i]) continue;
            const bam_pileup1_t *p = &plp[0][i];
            const bam1_t *b = p->b;
            int s = bio_strand(b);
            uint8_t nt = bam_seqi(bam_get_seq(b), p->qpos);
            int base_i = nt16_index(nt);
            /* Minus reads: stored SEQ is 5'->3' (SAM spec) but qpos walks
             * CIGAR order (left->right); complement the base into the
             * reference frame (parity with countmut 0.0.x + pysam pairs). */
            if (s == 1 && base_i < 4) base_i = 3 - base_i;
            /* router-assigned category slot (0..CM_CAT_MAX-1); g_a set in
             * the selection loop for every kept candidate. */
            int cat = w->g_a[i];
            if (cat < 0) cat = 0;                 /* defensive */
            site.cnt[s][cat][base_i]++;
        }

        /* ---------- emit ---------- */
        int emit_plus = cfg->strand_process != CM_STRAND_REVERSE;
        int emit_minus = cfg->strand_process != CM_STRAND_FORWARD;
        /* -p site filter (per-strand): keep the site if either strand passes,
         * then emit only the passing strand(s). */
        int smask = expr_pile_apply_strands(w->expr, cfg, w, hdr->target_name[tid],
                                            &site, pos, ref_ch);
        if (smask == 0) continue;
        emit_plus = emit_plus && (smask & 1);
        emit_minus = emit_minus && (smask & 2);
        emit_site(w, cfg, hdr, fp, tid, pos, ref_ch, &site, emit_plus, emit_minus);
    }
    free(n_plp); free(plp);
    bam_mplp_destroy(mplp);
    if (aux.itr) bam_itr_destroy(aux.itr);
}

/* ===========================================================================
 * Read-walk engine
 *
 * Walks the BAM read by read (like countmut's original core / the Python
 * engine_readwalk): each read's CIGAR is walked to the reference, matched
 * bases are deduplicated by (ref_pos, qname) with the (mapq, read1, qual)
 * preference tuple, and deletions / ref-skips / filter-failures are tallied
 * per read.  The result is flushed through the same emit_site() as the
 * pileup engine, so the two walks are byte-identical.
 * ======================================================================== */

/* dedup hash: (0-based pos, qname-id) -> index into the winner array.
 * qname ids come from a per-region qname->id table, so the (pos,qname) overlap
 * dedup is unchanged but keys are plain integers (no strdup per entry). */
typedef struct { int64_t pos; int qid; } posq_key;
/* SplitMix64-style finalizer: klib's kh_int64_hash_func has weak low bits, and
 * our keys (pos small, qid increasing in insertion order) made those low bits
 * near-monotonic, degrading klib's open addressing to O(n) per insert (O(n^2)
 * total) as the hash grew at deep sites.  This mix scrambles the low bits so
 * insertion order no longer correlates with bucket order. */
static inline khint_t posq_hash(posq_key k) {
    uint64_t z = (uint64_t)k.pos * 0x9E3779B97F4A7C15ull
               ^ (uint64_t)(k.qid + 1) * 0xBF58476D1CE4E5B9ull;
    z ^= z >> 30; z *= 0xBF58476D1CE4E5B9ull;
    z ^= z >> 27; z *= 0x94D049BB133111EBull;
    z ^= z >> 31;
    return (khint_t)z;
}
static inline int posq_equal(posq_key a, posq_key b) {
    return a.pos == b.pos && a.qid == b.qid;
}
KHASH_INIT(posq, posq_key, int, 1, posq_hash, posq_equal)
KHASH_INIT(qn2id, char *, int, 1, kh_str_hash_func, kh_str_hash_equal)

/* pos -> sitemap slot; strong mixer (see the posq_hash note -- klib's low bits
 * are weak and the direct path inserts positions in monotonic order) */
static inline khint_t posi_hash(khint64_t key) {
    uint64_t z = (uint64_t)key * 0x9E3779B97F4A7C15ull;
    z ^= z >> 30; z *= 0xBF58476D1CE4E5B9ull;
    z ^= z >> 27; z *= 0x94D049BB133111EBull;
    z ^= z >> 31;
    return (khint_t)z;
}
KHASH_INIT(posi, khint64_t, int, 1, posi_hash, kh_int64_hash_equal)

/* winner of a (pos,qname) dedup bucket */
typedef struct { int mapq, r1, qual, strand, base, slot; } rw_w;

/* growable map pos -> site_t (one entry per visited reference position).
 * When the fast-dna fast path is active, `arr` is a dense fixed array over
 * [arr_beg, arr_beg+arr_len) and sm_get() is an O(1) index instead of a hash
 * lookup; otherwise the hash map is used. */
typedef struct {
    khash_t(posi) *pm;
    int64_t *spos;
    site_t *st;
    int n, cap;
    site_t *arr; int64_t arr_beg; int arr_len;   /* fast-dna dense array */
} sitemap_t;

static void sm_init(sitemap_t *m, int fast_dna, int64_t beg, int64_t end) {
    m->pm = NULL; m->spos = NULL; m->st = NULL; m->n = m->cap = 0;
    m->arr = NULL; m->arr_beg = beg; m->arr_len = 0;
    if (fast_dna && end > beg) {
        m->arr_len = (int)(end - beg);
        m->arr = (site_t *)calloc((size_t)m->arr_len, sizeof(site_t));
        m->arr_beg = beg;
    } else {
        m->pm = kh_init(posi);
    }
}
static void sm_free(sitemap_t *m) {
    if (m->pm) kh_destroy(posi, m->pm);
    free(m->spos); free(m->st);
    if (m->arr) free(m->arr);
    memset(m, 0, sizeof(*m));
}
/* NOTE: the returned pointer is only valid until the next sm_get() (a later
 * insert may realloc), so it must be used immediately and never retained. */
static site_t *sm_get(sitemap_t *m, int64_t pos) {
    if (m->arr) {
        int64_t i = pos - m->arr_beg;
        if (i >= 0 && i < m->arr_len) return &m->arr[i];
        return NULL;   /* out of region: caller should not hit this */
    }
    khint_t k = kh_get(posi, m->pm, pos);
    if (k != kh_end(m->pm)) return &m->st[kh_val(m->pm, k)];
    int ret; k = kh_put(posi, m->pm, pos, &ret);
    if (m->n == m->cap) {
        m->cap = m->cap ? m->cap * 2 : 64;
        m->st = (site_t *)realloc(m->st, (size_t)m->cap * sizeof(site_t));
        m->spos = (int64_t *)realloc(m->spos, (size_t)m->cap * sizeof(int64_t));
    }
    int idx = m->n++;
    kh_val(m->pm, k) = idx;
    m->spos[idx] = pos;
    site_zero(&m->st[idx]);
    return &m->st[idx];
}

typedef struct { int64_t pos; int idx; } site_ord;
static int cmp_site_ord(const void *a, const void *b) {
    int64_t x = ((const site_ord *)a)->pos, y = ((const site_ord *)b)->pos;
    return (x > y) - (x < y);
}

/* ---- read-walk fast-path helpers -------------------------------------------
 * The slow path walks every aligned base; for indel-free reads the fast path
 * jumps straight to the (sorted) target positions.  Both funnel into
 * rw_add_base(), so dedup/quality decisions are byte-identical. */

/* Insert/update one (pos,qid) winner in the posq dedup hash.  Shared by
 * rw_add_base() and the duplicate-qname promotion path (frag_ovl_promote). */
static rw_w *posq_add(khash_t(posq) *h, rw_w *wins, int *wins_cap, int *wins_n,
                      int64_t ref_pos, int qid, int mapq, int r1, int qual,
                      int strand, int base, int slot) {
    posq_key key = { ref_pos, qid };
    int r;
    khint_t kh = kh_put(posq, h, key, &r);
    if (r) {   /* new (pos,qname) */
        if (*wins_n == *wins_cap) {
            *wins_cap = *wins_cap ? *wins_cap * 2 : 64;
            wins = (rw_w *)realloc(wins, (size_t)*wins_cap * sizeof(rw_w));
        }
        int idx = (*wins_n)++;
        wins[idx].mapq = mapq; wins[idx].r1 = r1; wins[idx].qual = qual;
        wins[idx].strand = strand; wins[idx].base = base; wins[idx].slot = slot;
        kh_val(h, kh) = idx;
    } else {
        int j = kh_val(h, kh);
        if (better(mapq, r1, qual, wins[j].mapq, wins[j].r1, wins[j].qual)) {
            wins[j].mapq = mapq; wins[j].r1 = r1; wins[j].qual = qual;
            wins[j].strand = strand; wins[j].base = base; wins[j].slot = slot;
        }
    }
    return wins;
}

/* Add one matched base to the (pos,qname) dedup table.  Applies the mutation
 * target gate (skipped when `already_target`), trim (is_internal), the -e
 * filter and the (mapq,read1,qual) preference.  When `direct` is set the base
 * sits outside any mate-overlap region (or the read is single-end / only one
 * mate covers), so it is counted straight into the site with no hash -- the
 * fragment contributes that base exactly once either way, keeping output
 * identical for the intended (overlapping-mate) dedup.  `qid` is the read's
 * qname id (resolved once per read).  Returns the (possibly reallocated) wins. */
static rw_w *rw_add_base(worker_t *w, const cm_config *cfg, bam_hdr_t *hdr, int tid,
                         const bam1_t *b, int s, int64_t ref_pos, uint32_t qpos,
                         int qid, int already_target, int direct, sitemap_t *sm,
                         khash_t(posq) *h, rw_w *wins, int *wins_cap, int *wins_n) {
    uint32_t qlen = b->core.l_qseq;
    if (qpos >= qlen) return wins;
    /* -e router: category slot (-1 = drop).  Read-constant -e is applied once
     * per read by the caller (see count_interval_readwalk); here we route the
     * per-base (non-read-constant) case only, which is what a bq/qpos-based
     * router needs.  Default slot is 0. */
    int rslot = 0;
    if (w->expr && cm_expr_has_read(w->expr) && !cm_expr_read_constant(w->expr)) {
        rslot = cm_expr_route(w->expr, b, hdr->target_name[tid],
                              (b->core.mtid >= 0 && b->core.mtid < hdr->n_targets)
                                  ? hdr->target_name[b->core.mtid] : "",
                              (int)qpos, s ? -1 : 1,
                              (ref_pos >= 0 && ref_pos < w->chr_len) ? w->chr_seq[ref_pos] : 'N');
        if (rslot < 0) return wins;
    }
    uint8_t nt = bam_seqi(bam_get_seq(b), qpos);
    int base_i = nt16_index(nt);
    /* Minus reads: stored SEQ is 5'->3' (SAM spec) but qpos walks CIGAR order
     * (left->right); complement the base into the reference frame (parity
     * with countmut 0.0.x + pysam pairs). */
    if (s == 1 && base_i < 4) base_i = 3 - base_i;
    int qual = (int)bam_get_qual(b)[qpos];
    if (direct) {
        int cat = rslot;
        sm_get(sm, ref_pos)->cnt[s][cat][base_i]++;
        return wins;
    }
    int mapq = (int)b->core.qual;
    int r1 = (b->core.flag & BAM_FREAD1) ? 1 : 0;
    return posq_add(h, wins, wins_cap, wins_n, ref_pos, qid, mapq, r1, qual, s, base_i, rslot);
}

/* Per-fragment overlap buffer for fast-dna mode.
 *
 * In fast-dna mode the global (pos,qname) dedup hash becomes the bottleneck
 * at deep coverage (tens to hundreds of millions of inserts, dominated by
 * cache misses).  For paired reads whose mates overlap we instead buffer the
 * overlap bases per fragment and resolve the winner once both mates are seen.
 * The overlap reference interval is contiguous, so we store winners by offset
 * inside that interval -- no hash table on the per-base hot path.
 *
 * For a fragment with read1 at pos1 and read2 at pos2, the outer template
 * length is isize = pos2 + qlen - pos1.  When isize < 2*qlen the mates overlap
 * on reference [pos2, pos1+qlen), length L = 2*qlen - isize.  Both mates map
 * their overlap bases to the same L reference positions, so a per-offset
 * comparison picks the (mapq,read1,qual) winner exactly like the posq path. */
typedef struct {
    int seen;              /* bitmask: 1=read1 seen, 2=read2 seen */
    int n_read;            /* reads seen with this qname (detect >2 duplicates) */
    int L;                 /* overlap length on reference */
    int64_t ref_start;     /* 0-based start of overlap on reference */
    uint8_t *has[2];       /* has[0]=read1, has[1]=read2 */
    rw_w *win[2];          /* winner candidate per offset per mate */
} frag_ovl_t;

static frag_ovl_t *frag_ovl_new(int L, int64_t ref_start) {
    frag_ovl_t *f = (frag_ovl_t *)calloc(1, sizeof(frag_ovl_t));
    if (L < 1) L = 1;
    f->L = L; f->ref_start = ref_start; f->seen = 0;
    for (int m = 0; m < 2; ++m) {
        f->has[m] = (uint8_t *)calloc((size_t)L, sizeof(uint8_t));
        f->win[m] = (rw_w *)malloc((size_t)L * sizeof(rw_w));
    }
    return f;
}
static void frag_ovl_free(frag_ovl_t *f) {
    if (!f) return;
    for (int m = 0; m < 2; ++m) { free(f->has[m]); free(f->win[m]); }
    free(f);
}
/* Ensure the buffer covers [ref_start, ref_start+L).  Expand and shift if the
 * second mate reports a slightly different interval (can happen with clipping). */
static void frag_ovl_ensure(frag_ovl_t *f, int L2, int64_t ref_start2) {
    int64_t end1 = f->ref_start + f->L;
    int64_t end2 = ref_start2 + L2;
    int64_t new_start = f->ref_start < ref_start2 ? f->ref_start : ref_start2;
    int64_t new_end = end1 > end2 ? end1 : end2;
    int new_L = (int)(new_end - new_start);
    if (new_L <= f->L && new_start == f->ref_start) return;
    for (int m = 0; m < 2; ++m) {
        uint8_t *nh = (uint8_t *)calloc((size_t)new_L, sizeof(uint8_t));
        rw_w *nw = (rw_w *)malloc((size_t)new_L * sizeof(rw_w));
        if (f->has[m]) {
            int off = (int)(f->ref_start - new_start);
            memcpy(nh + off, f->has[m], (size_t)f->L * sizeof(uint8_t));
            memcpy(nw + off, f->win[m], (size_t)f->L * sizeof(rw_w));
            free(f->has[m]); free(f->win[m]);
        }
        f->has[m] = nh; f->win[m] = nw;
    }
    f->L = new_L; f->ref_start = new_start;
}
static void frag_ovl_store(frag_ovl_t *f, int mate, int64_t ref_pos,
                           int mapq, int r1, int qual, int strand, int base, int slot) {
    int off = (int)(ref_pos - f->ref_start);
    if (off < 0 || off >= f->L) {
        /* Clipping/indel can shift the actual aligned overlap relative to the
         * geometry-derived interval; expand (and shift) the buffer so every
         * overlap base lands inside it. */
        frag_ovl_ensure(f, 1, ref_pos);
        off = (int)(ref_pos - f->ref_start);
        if (off < 0 || off >= f->L) return;
    }
    f->has[mate][off] = 1;
    f->win[mate][off].mapq = mapq;
    f->win[mate][off].r1 = r1;
    f->win[mate][off].qual = qual;
    f->win[mate][off].strand = strand;
    f->win[mate][off].base = base;
    f->win[mate][off].slot = slot;
}
static void frag_ovl_flush(frag_ovl_t *f, sitemap_t *sm, int64_t beg, int64_t end) {
    for (int off = 0; off < f->L; ++off) {
        int64_t pos = f->ref_start + off;
        if (pos < beg || pos >= end) continue;
        rw_w *w = NULL;
        if (f->has[0][off] && f->has[1][off]) {
            rw_w *a = &f->win[0][off], *b = &f->win[1][off];
            w = better(a->mapq, a->r1, a->qual, b->mapq, b->r1, b->qual) ? a : b;
        } else if (f->has[0][off]) {
            w = &f->win[0][off];
        } else if (f->has[1][off]) {
            w = &f->win[1][off];
        }
        if (w) {
            site_t *st = sm_get(sm, pos);
            if (st) st->cnt[w->strand][w->slot][w->base]++;
        }
    }
}
KHASH_INIT(fovl, int, frag_ovl_t *, 1, kh_int_hash_func, kh_int_hash_equal)
KHASH_INIT(dupeset, int, int, 0, kh_int_hash_func, kh_int_hash_equal)

/* A qname appearing >2 times (split/duplicate alignments) cannot be handled by
 * the two-slot per-fragment buffer; migrate its buffered bases into the global
 * posq hash so all reads with that qname dedup against each other exactly as
 * the original posq path would.  Returns the (possibly reallocated) wins. */
static rw_w *frag_ovl_promote(frag_ovl_t *f, int qid, khash_t(posq) *h,
                              rw_w *wins, int *wins_cap, int *wins_n) {
    for (int mate = 0; mate < 2; ++mate) {
        for (int off = 0; off < f->L; ++off) {
            if (!f->has[mate][off]) continue;
            int64_t ref_pos = f->ref_start + off;
            rw_w *w = &f->win[mate][off];
            wins = posq_add(h, wins, wins_cap, wins_n, ref_pos, qid,
                            w->mapq, w->r1, w->qual, w->strand, w->base, w->slot);
        }
    }
    return wins;
}

static void count_interval_readwalk(worker_t *w, const cm_config *cfg, bam_hdr_t *hdr,
                                    FILE *fp, int tid, int beg, int end) {
    aux_t aux;
    aux.fp = w->fp; aux.hdr = hdr;
    aux.beg = beg; aux.end = end;
    aux.itr = w->idx ? sam_itr_queryi(w->idx, tid, beg, end) : NULL;

    if (w->last_tid != tid) load_chr_seq(w, hdr, tid);

    khash_t(posq) *h = kh_init(posq);
    khash_t(qn2id) *qnids = kh_init(qn2id);
    khash_t(fovl) *fovl_h = kh_init(fovl);
    khash_t(dupeset) *dupe_qids = kh_init(dupeset);
    int qname_n = 0;
    rw_w *wins = NULL; int wins_cap = 0, wins_n = 0;
    /* Dense fixed-array sitemap (O(1) index) is the default for strandless
     * counting.  It supports the -e read filter: a simple "bq >= N" filter is
     * applied as a fast integer check, any other -e expression is evaluated
     * per base via the Lua router.  The general hash-map sitemap is used only
     * when the counting needs motif/BED/target/indel features. */
    int fast_dna = cfg->strandless &&
                   !cfg->pile_expr && cfg->pad == 0 && !cfg->bedfile &&
                   !cfg->exclude && cfg->target_base < 0 && !cfg->count_indels;
    sitemap_t sm; sm_init(&sm, fast_dna, beg, end);
    if (cfg->verbose && fast_dna) fprintf(stderr, "[countmut] fast-dna region %d:%d-%d\n", tid, beg, end);

    bam1_t *b = bam_init1();
    int ret;
    while ((ret = (aux.itr ? sam_itr_next(aux.fp, aux.itr, b) : sam_read1(aux.fp, aux.hdr, b))) >= 0) {
        int s = bio_strand(b);
        if (cfg->strand_process == CM_STRAND_FORWARD && s != 0) continue;
        if (cfg->strand_process == CM_STRAND_REVERSE && s != 1) continue;
        if (read_fails(cfg, b)) {
            /* filter-failure: tally `fail` at every covered reference position */
            const uint32_t *cig = bam_get_cigar(b);
            int64_t rcur = b->core.pos;
            for (int i = 0; i < b->core.n_cigar; ++i) {
                int op = (int)bam_cigar_op_p(&cig[i]); int len = (int)bam_cigar_oplen_p(&cig[i]);
                if (op == 0 || op == 2 || op == 3 || op == 7 || op == 8) { /* consumes ref */
                    if (rcur < end && rcur + len > beg) {
                        int64_t lo = rcur > beg ? rcur : beg;
                        int64_t hi = rcur + len < end ? rcur + len : end;
                        for (int64_t p = lo; p < hi; ++p)
                            sm_get(&sm, p)->fail[s]++;
                    }
                    rcur += len;
                }
            }
            continue;
        }

        uint32_t qlen = b->core.l_qseq;
        if (qlen == 0) continue;
        const uint32_t *cig = bam_get_cigar(b);

        /* Solo-vs-overlap for the overlap dedup: a base is "direct" (counted
         * straight into the site) unless this read's mate can also cover it.
         * Only single-end reads / mates on another contig are guaranteed solo
         * (all-direct).  Overlapping mates with known geometry use the hybrid
         * path; anything uncertain (mpos<0 / TLEN unusable) keeps the exact
         * hash for every base.
         *   read1 overlap: qpos >= ins - r     read2: qpos < 2r - ins */
        int ovl = -1, olo = 0, ohi = (int)qlen;
        if (!(b->core.flag & BAM_FPAIRED)) ovl = 0;                       /* single-end: no mate */
        else if (b->core.mtid >= 0 && b->core.mtid != b->core.tid) ovl = 0; /* mate known elsewhere */
        else if (b->core.mpos >= 0 && b->core.mtid == b->core.tid) {
            int mdiff = abs((int)b->core.mpos - (int)b->core.pos);
            if (mdiff >= (int)qlen) ovl = 0;   /* mate spans are disjoint (|mpos-pos| >= read length):
                                               * no reference position has two covers -> count direct,
                                               * skipping the (pos,qname) dedup hash entirely */
            else {
                /* mates can overlap (starts < one read length apart): dedup the
                 * overlap window via the hash; known fragment geometry gives the
                 * exact window (hybrid), otherwise use the exact hash path */
                int insv = (int)b->core.isize; if (insv < 0) insv = -insv;
                if (insv > 0 && insv < 2 * (int)qlen) {
                    ovl = 1;
                    if (b->core.flag & BAM_FREAD1) olo = insv - (int)qlen;
                    else ohi = 2 * (int)qlen - insv;
                }
                /* else ovl stays -1 -> exact hash (safe fallback) */
            }
        }
        /* else (paired but mate position unknown: mpos<0/mtid<0) ovl stays -1
         * -> every base goes through the hash = exact (never direct) */
        /* fast-dna overlap geometry (pre-compute once per read). */
        int L_ovl = 0;
        int64_t ref_start_ovl = 0;
        if (fast_dna && ovl == 1) {
            if (b->core.flag & BAM_FREAD1) {
                L_ovl = (int)qlen - olo;
                ref_start_ovl = (int64_t)b->core.pos + olo;
            } else {
                L_ovl = ohi;
                ref_start_ovl = (int64_t)b->core.pos;
            }
            if (L_ovl < 1) L_ovl = 1;
        }
        /* qname -> id: only needed for the (pos,qname) overlap dedup hash.
         * In the fast-dna path a non-overlapping read (ovl==0) is counted
         * directly into the dense array, so we skip the strdup/hash entirely. */
        int qid = 0;
        if (!(fast_dna && ovl == 0)) {
            const char *qnbuf = bam_get_qname(b);
            khint_t qk = kh_get(qn2id, qnids, qnbuf);
            if (qk == kh_end(qnids)) {
                int r; char *cp = strdup(qnbuf); qk = kh_put(qn2id, qnids, cp, &r);
                qid = qname_n++; kh_val(qnids, qk) = qid;
            } else {
                qid = kh_val(qnids, qk);
            }
        }
        /* read-constant -e filter: evaluate ONCE per read (not per base) and
         * skip the whole read when it fails.  rw_add_base() then skips the
         * per-base -e call for these (its decision is already known). */
        if (w->expr && cm_expr_has_read(w->expr) && cm_expr_read_constant(w->expr)) {
            const char *mrn = (b->core.mtid >= 0 && b->core.mtid < hdr->n_targets)
                                  ? hdr->target_name[b->core.mtid] : "";
            if (!cm_expr_read(w->expr, b, hdr->target_name[tid], mrn, 0, s ? -1 : 1, 'N'))
                continue;
        }
        /* direct(ref_pos,qpos) = counts straight into the site (no dedup hash) */
#define RW_DIRECT(_qpos) ((ovl == 0) || (ovl == 1 && ((int)(_qpos) < olo || (int)(_qpos) >= ohi)))

        {
            uint32_t qcur = 0;
            int64_t rcur = b->core.pos;
            for (int i = 0; i < b->core.n_cigar; ++i) {
                int op = (int)bam_cigar_op_p(&cig[i]); int len = (int)bam_cigar_oplen_p(&cig[i]);
                switch (op) {
                case 0: case 7: case 8: /* M, =, X -- matched bases */
                    if (fast_dna && ovl == 0) {
                        /* fast-dna inline path: ovl==0 so every base is direct,
                         * count straight into the dense array (O(1) index), no
                         * rw_add_base call / no (pos,qname) dedup hash. */
                        for (int k = 0; k < len; ++k) {
                            int64_t ref_pos = rcur + k;
                            if (ref_pos < beg || ref_pos >= end) continue;
                            uint32_t qpos = qcur + (uint32_t)k;
                            if (qpos >= qlen) break;
                            uint8_t nt = bam_seqi(bam_get_seq(b), qpos);
                            int base_i = nt16_index(nt);
                            if (s == 1 && base_i < 4) base_i = 3 - base_i;
                            int cat = 0;
                            if (w->expr && cm_expr_has_read(w->expr) && !cm_expr_read_constant(w->expr)) {
                                cat = cm_expr_route(w->expr, b, hdr->target_name[tid],
                                    (b->core.mtid >= 0 && b->core.mtid < hdr->n_targets) ? hdr->target_name[b->core.mtid] : "",
                                    (int)qpos, s ? -1 : 1,
                                    (ref_pos >= 0 && ref_pos < w->chr_len) ? w->chr_seq[ref_pos] : 'N');
                                if (cat < 0) continue;
                            }
                            site_t *st = sm_get(&sm, ref_pos);
                            if (st) st->cnt[s][cat][base_i]++;
                        }
                    } else if (fast_dna && ovl == 1) {
                        /* fast-dna overlap path: buffer per-fragment, no posq hash. */
                        int mate = (b->core.flag & BAM_FREAD2) ? 1 : 0;  /* 0=R1, 1=R2 */
                        int mate_bit = 1 << mate;
                        int mapq = (int)b->core.qual;
                        int r1f = (b->core.flag & BAM_FREAD1) ? 1 : 0;
                        khint_t dk = kh_get(dupeset, dupe_qids, qid);
                        if (dk != kh_end(dupe_qids)) {
                            /* duplicate qname: route through posq so all reads
                             * with this qname dedup against each other. */
                            for (int k = 0; k < len; ++k) {
                                int64_t ref_pos = rcur + k;
                                if (ref_pos < beg || ref_pos >= end) continue;
                                uint32_t qpos = qcur + (uint32_t)k;
                                if (qpos >= qlen) break;
                                wins = rw_add_base(w, cfg, hdr, tid, b, s, ref_pos, qpos,
                                                   qid, 0, RW_DIRECT(qpos), &sm,
                                                   h, wins, &wins_cap, &wins_n);
                            }
                        } else {
                            khint_t fk = kh_get(fovl, fovl_h, qid);
                            frag_ovl_t *fo;
                            int fnew;
                            if (fk == kh_end(fovl_h)) {
                                fk = kh_put(fovl, fovl_h, qid, &fnew);
                                fo = frag_ovl_new(L_ovl, ref_start_ovl);
                                fo->n_read = 1;
                                kh_val(fovl_h, fk) = fo;
                            } else {
                                fo = kh_val(fovl_h, fk);
                                fo->n_read++;
                            }
                            if (fo->n_read > 2) {
                                /* third+ read with this qname: promote to posq. */
                                wins = frag_ovl_promote(fo, qid, h, wins, &wins_cap, &wins_n);
                                kh_del(fovl, fovl_h, fk);
                                frag_ovl_free(fo);
                                int r; kh_put(dupeset, dupe_qids, qid, &r);
                                for (int k = 0; k < len; ++k) {
                                    int64_t ref_pos = rcur + k;
                                    if (ref_pos < beg || ref_pos >= end) continue;
                                    uint32_t qpos = qcur + (uint32_t)k;
                                    if (qpos >= qlen) break;
                                    wins = rw_add_base(w, cfg, hdr, tid, b, s, ref_pos, qpos,
                                                       qid, 0, RW_DIRECT(qpos), &sm,
                                                       h, wins, &wins_cap, &wins_n);
                                }
                            } else {
                                frag_ovl_ensure(fo, L_ovl, ref_start_ovl);
                                fo->seen |= mate_bit;
                                for (int k = 0; k < len; ++k) {
                                    int64_t ref_pos = rcur + k;
                                    if (ref_pos < beg || ref_pos >= end) continue;
                                    uint32_t qpos = qcur + (uint32_t)k;
                                    if (qpos >= qlen) break;
                                    int in_ovl = ((int)qpos >= olo && (int)qpos < ohi);
                                    uint8_t nt = bam_seqi(bam_get_seq(b), qpos);
                                    int base_i = nt16_index(nt);
                                    if (s == 1 && base_i < 4) base_i = 3 - base_i;
                                    int cat = 0;
                                    if (w->expr && cm_expr_has_read(w->expr) && !cm_expr_read_constant(w->expr)) {
                                        cat = cm_expr_route(w->expr, b, hdr->target_name[tid],
                                            (b->core.mtid >= 0 && b->core.mtid < hdr->n_targets) ? hdr->target_name[b->core.mtid] : "",
                                            (int)qpos, s ? -1 : 1,
                                            (ref_pos >= 0 && ref_pos < w->chr_len) ? w->chr_seq[ref_pos] : 'N');
                                        if (cat < 0) continue;
                                    }
                                    if (!in_ovl) {
                                        site_t *st = sm_get(&sm, ref_pos);
                                        if (st) st->cnt[s][cat][base_i]++;
                                    } else {
                                        int qual = (int)bam_get_qual(b)[qpos];
                                        frag_ovl_store(fo, mate, ref_pos, mapq, r1f, qual, s, base_i, cat);
                                    }
                                }
                                if (fo->seen == 3) {
                                    frag_ovl_flush(fo, &sm, beg, end);
                                    kh_del(fovl, fovl_h, fk);
                                    frag_ovl_free(fo);
                                }
                            }
                        }
                    } else {
                        for (int k = 0; k < len; ++k) {
                            int64_t ref_pos = rcur + k;
                            if (ref_pos < beg || ref_pos >= end) continue;
                            uint32_t qpos = qcur + (uint32_t)k;
                            if (qpos >= qlen) break;
                            wins = rw_add_base(w, cfg, hdr, tid, b, s, ref_pos, qpos,
                                               qid, 0, RW_DIRECT(qpos), &sm,
                                               h, wins, &wins_cap, &wins_n);
                        }
                    }
                    qcur += (uint32_t)len; rcur += len;
                    break;
                case 1: case 4: qcur += (uint32_t)len; break; /* I, S consume query */
                case 2: case 3: { /* D, N -- deletion / ref-skip */
                    int is_del = (op == 2);
                    if (rcur < end && rcur + len > beg) {
                        int64_t lo2 = rcur > beg ? rcur : beg;
                        int64_t hi2 = rcur + len < end ? rcur + len : end;
                        for (int64_t p = lo2; p < hi2; ++p)
                            if (is_del) sm_get(&sm, p)->del[s]++; else sm_get(&sm, p)->refskip[s]++;
                    }
                    rcur += len;
                    break;
                }
                default: break; /* H, P consume nothing */
                }
            }
        }
    }

    /* flush dedup winners into the per-position sites */
    for (khint_t k = kh_begin(h); k != kh_end(h); ++k) {
        if (!kh_exist(h, k)) continue;
        int64_t pos = kh_key(h, k).pos;
        const rw_w *win = &wins[kh_val(h, k)];
        int cat = win->slot;
        site_t *st = sm_get(&sm, pos);
        st->cnt[win->strand][cat][win->base]++;
    }

    /* flush any orphan overlap fragments (mate missing / outside region) */
    for (khint_t k = kh_begin(fovl_h); k != kh_end(fovl_h); ++k) {
        if (!kh_exist(fovl_h, k)) continue;
        frag_ovl_t *fo = kh_val(fovl_h, k);
        frag_ovl_flush(fo, &sm, beg, end);
        frag_ovl_free(fo);
    }

    /* emit sites in position order */
    const int emit_plus = cfg->strand_process != CM_STRAND_REVERSE;
    const int emit_minus = cfg->strand_process != CM_STRAND_FORWARD;
    if (sm.arr) {
        /* fast-dna dense array: positions are contiguous, no sort needed. */
        for (int i = 0; i < sm.arr_len; ++i) {
            int64_t pos = sm.arr_beg + i;
            if (pos < 0 || pos >= w->chr_len) continue;
            char ref_ch = w->chr_seq[pos];
            if (w->inc_bed && !bed_overlap(w->inc_bed, hdr->target_name[tid], (int)pos, (int)pos + 1)) continue;
            if (w->exc_bed && bed_overlap(w->exc_bed, hdr->target_name[tid], (int)pos, (int)pos + 1)) continue;
            int smask = expr_pile_apply_strands(w->expr, cfg, w, hdr->target_name[tid],
                                                &sm.arr[i], pos, ref_ch);
            if (smask == 0) continue;
            emit_site(w, cfg, hdr, fp, tid, (int)pos, ref_ch, &sm.arr[i],
                      emit_plus && (smask & 1), emit_minus && (smask & 2));
        }
    } else {
        site_ord *ord = (site_ord *)malloc((size_t)(sm.n ? sm.n : 1) * sizeof(*ord));
        for (int i = 0; i < sm.n; ++i) { ord[i].pos = sm.spos[i]; ord[i].idx = i; }
        qsort(ord, (size_t)sm.n, sizeof(*ord), cmp_site_ord);
        for (int i = 0; i < sm.n; ++i) {
            int64_t pos = ord[i].pos;
            if (pos < 0 || pos >= w->chr_len) continue;
            char ref_ch = w->chr_seq[pos];
            if (w->inc_bed && !bed_overlap(w->inc_bed, hdr->target_name[tid], (int)pos, (int)pos + 1)) continue;
            if (w->exc_bed && bed_overlap(w->exc_bed, hdr->target_name[tid], (int)pos, (int)pos + 1)) continue;
            int smask = expr_pile_apply_strands(w->expr, cfg, w, hdr->target_name[tid],
                                                &sm.st[ord[i].idx], pos, ref_ch);
            if (smask == 0) continue;
            emit_site(w, cfg, hdr, fp, tid, (int)pos, ref_ch, &sm.st[ord[i].idx],
                      emit_plus && (smask & 1), emit_minus && (smask & 2));
        }
        free(ord);
    }

    /* cleanup */
    for (khint_t k = kh_begin(qnids); k != kh_end(qnids); ++k)
        if (kh_exist(qnids, k)) free((void *)(uintptr_t)kh_key(qnids, k));
    kh_destroy(qn2id, qnids);
    kh_destroy(posq, h);
    kh_destroy(fovl, fovl_h);
    kh_destroy(dupeset, dupe_qids);
    free(wins);
    sm_free(&sm);
    bam_destroy1(b);
    if (aux.itr) bam_itr_destroy(aux.itr);
}

typedef struct { work_t *s; int wi; } targ_t;

/* current resident memory (MB) of this process, for --verbose progress */
static long cur_rss_mb(void) {
    long pages = 0, size = 0;
    FILE *f = fopen("/proc/self/statm", "r");
    if (f) { if (fscanf(f, "%ld %ld", &size, &pages) != 2) pages = 0; fclose(f); }
    return pages * ((long)sysconf(_SC_PAGESIZE) / 1024) / 1024;
}

static void *thread_main(void *arg) {
    targ_t *ta = (targ_t *)arg;
    worker_t *w = &ta->s->workers[ta->wi];
    work_t *s = ta->s;
    FILE *wf = s->files[ta->wi];   /* one temp file per worker (bounded fd count) */
    const int step = s->nregions > 100 ? s->nregions / 100 : 1;
    for (;;) {
        int i = __sync_fetch_and_add(&s->next, 1);   /* dynamic claim, keeps deep
                                                      * bins from serializing on
                                                      * one static slice */
        if (i >= s->nregions) break;
        long off = ftell(wf);
        if (s->cfg->engine == CM_ENGINE_READWALK)
            count_interval_readwalk(w, s->cfg, s->hdr, wf,
                                    s->regions[i].tid, s->regions[i].beg, s->regions[i].end);
        else
            count_interval(w, s->cfg, s->hdr, wf,
                           s->regions[i].tid, s->regions[i].beg, s->regions[i].end);
        s->spans[i].worker = ta->wi;
        s->spans[i].off = off;
        s->spans[i].len = ftell(wf) - off;
        if (s->cfg->verbose) {
            int done = __sync_add_and_fetch(&s->done, 1);
            if (done % step == 0 || done == s->nregions)
                fprintf(stderr, "[countmut] %d/%d regions (%.1f%%) done  rss=%ldMB\n",
                        done, s->nregions, 100.0 * done / s->nregions, cur_rss_mb());
        }
    }
    return NULL;
}

/* Build genomic bins so there are roughly 4*tasks bins per genome. */
static region_t *build_regions(bam_hdr_t *hdr, int threads, const char *region, int *n_out) {
    region_t *regs = NULL; int n = 0, cap = 0;
    if (region) {
        char chr[256]; long st = 0, en = 0; char *s;
        strncpy(chr, region, 255); chr[255] = 0;
        s = strchr(chr, ':');
        if (!s) return NULL;
        *s = 0; sscanf(s + 1, "%ld-%ld", &st, &en);
        /* Samtools-style 1-based inclusive region: st -> 0-based start (st-1),
         * en stays as the 0-based exclusive end.  Matches the Python wrapper. */
        if (st > 0) --st;
        if (st < 0) st = 0;
        int tid = bam_name2id(hdr, chr);
        if (tid < 0) return NULL;
        regs = (region_t *)malloc(sizeof(region_t)); regs[0].tid = tid; regs[0].beg = (int)st; regs[0].end = (int)en;
        *n_out = 1; return regs;
    }
    long total = 0; int nt = hdr->n_targets;
    for (int i = 0; i < nt; ++i) total += hdr->target_len[i];
    long bin_size = total / (4L * (threads < 1 ? 1 : threads));
    if (bin_size < 1) bin_size = 1;
    for (int i = 0; i < nt; ++i) {
        long len = hdr->target_len[i];
        for (long b = 0; b < len; b += bin_size) {
            if (n == cap) { cap = cap ? cap * 2 : 64; regs = (region_t *)realloc(regs, cap * sizeof(region_t)); }
            regs[n].tid = i; regs[n].beg = (int)b; regs[n].end = (int)(b + bin_size < len ? b + bin_size : len);
            ++n;
        }
    }
    *n_out = n; return regs;
}

/* ---- input format support -----------------------------------------------
 * BAM: native (BGZF + BAI index).
 * SAM (plain or gzipped): transcoded once to a temp BAM + BAI so the whole
 * indexed, multi-threaded pipeline is reused unchanged (identical output to
 * running on the equivalent BAM).
 * CRAM: NOT supported in this self-contained core (no CRAM codec); we fail
 * with a conversion hint instead of a confusing crash.
 * Returns 0=BAM, 1=SAM(transcoded), 2=CRAM(unsupported), -1=error. */
/* Input format detection.  Returns 0 = BAM/CRAM (handled natively by the
 * full htslib sam_open), 1 = SAM text (transcoded to a temp BAM), -1 = error.
 * CRAM is read natively: the core links the full htslib with a CRAM codec. */
static int detect_input_format(const char *path, int *is_sam) {
    unsigned char magic[4] = {0};
    gzFile gz = gzopen(path, "rb");
    if (gz == NULL) return -1;
    int n = (int)gzread(gz, magic, 4);
    gzclose(gz);
    if (n >= 4 && memcmp(magic, "BAM\1", 4) == 0) { *is_sam = 0; return 0; }
    if (n >= 4 && memcmp(magic, "CRAM", 4) == 0)   { *is_sam = 0; return 0; }
    *is_sam = 1; return 1;   /* SAM text (or empty -> header parse fails later with a clear error) */
}

static int transcode_sam_to_bam(const char *sam, char *tmp_bam, size_t cap) {
    samFile *in = sam_open(sam, "r");
    if (in == NULL) {
        fprintf(stderr, "[countmut] error: cannot open SAM input '%s'\n", sam);
        if (in) hts_close(in);
        return -1;
    }
    bam_hdr_t *hdr = sam_hdr_read(in);
    if (hdr == NULL) {
        fprintf(stderr, "[countmut] error: cannot parse SAM header from '%s'\n", sam);
        hts_close(in);
        return -1;
    }
    const char *td = getenv("TMPDIR");
    if (td == NULL || *td == '\0') td = "/tmp";
    char tpl[1024];
    snprintf(tpl, sizeof(tpl), "%s/countmut_sam_XXXXXX", td);
    int fd = mkstemp(tpl);
    if (fd < 0) {
        fprintf(stderr, "[countmut] error: cannot create temp file for SAM input\n");
        bam_hdr_destroy(hdr); hts_close(in);
        return -1;
    }
    close(fd);
    unlink(tpl);                          /* we only wanted the unique name */
    snprintf(tmp_bam, cap, "%s.bam", tpl);

    samFile *out = sam_open(tmp_bam, "wb");   /* BGZF-compressed BAM (indexable) */
    if (out == NULL) {
        fprintf(stderr, "[countmut] error: cannot write temp BAM '%s'\n", tmp_bam);
        bam_hdr_destroy(hdr); hts_close(in);
        unlink(tmp_bam);
        return -1;
    }
    sam_hdr_write(out, hdr);
    bam1_t *b = bam_init1();
    int nrec = 0;
    while (sam_read1(in, hdr, b) >= 0) {
        sam_write1(out, hdr, b);
        ++nrec;
    }
    sam_close(out);
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(in);
    if (bam_index_build(tmp_bam, 0) != 0) {
        fprintf(stderr, "[countmut] error: cannot index temp BAM '%s'\n", tmp_bam);
        return -1;
    }
    fprintf(stderr, "[countmut] input is SAM: converted %d records -> %s\n", nrec, tmp_bam);
    return 0;
}

int cm_run(const cm_config *cfg, const char *bam, const char *fa, const char *out_path, const char *region) {
    FILE *fp = (out_path && strcmp(out_path, "-") != 0) ? fopen(out_path, "w") : stdout;
    if (!fp) return 1;
    int is_sam = 0;
    char sam_tmp[1100] = {0};
    if (detect_input_format(bam, &is_sam) == 1) {
        /* SAM text: transcode to a temp BAM (the core reads BAM/CRAM natively). */
        if (transcode_sam_to_bam(bam, sam_tmp, sizeof(sam_tmp)) != 0) {
            if (fp != stdout) fclose(fp);
            return 3;
        }
        bam = sam_tmp;   /* the rest of the run operates on the temp BAM */
    }
    samFile *hfp = sam_open(bam, "r");
    if (!hfp) {
        fprintf(stderr, "[countmut] error: cannot open input file '%s'\n", bam);
        if (fp != stdout) fclose(fp);
        return 3;
    }
    hts_set_fai_filename(hfp, fa);   /* CRAM decode needs the reference */
    bam_hdr_t *hdr = sam_hdr_read(hfp);
    sam_close(hfp);
    if (!hdr) {
        fprintf(stderr, "[countmut] error: cannot read BAM header from '%s'\n", bam);
        if (fp != stdout) fclose(fp);
        return 3;
    }

    int nthreads = cfg->threads < 1 ? 1 : cfg->threads;
    /* Shared CRAM decode pool: all workers decode CRAM through one pool so
     * the per-file decoder threads never contend with the worker threads.
     * BAM reads are single-threaded anyway, so this only helps CRAM. */
    htsThreadPool tpool = {NULL, 0};
    int use_pool = 0;
    if (nthreads > 1) {
        tpool.pool = hts_tpool_init(nthreads);
        if (tpool.pool) { tpool.qsize = nthreads * 2; use_pool = 1; }
    }
    int nregions = 0;
    region_t *regs = build_regions(hdr, nthreads, region, &nregions);
    if (!regs) { bam_hdr_destroy(hdr); if (fp != stdout) fclose(fp); return 4; }
    if (nregions > 1 && nthreads > nregions) nthreads = nregions;

    write_header(fp, cfg);

    /* bounded temp files: one per worker (contiguous region slices), not one
     * per region -- an all-contigs run can have tens of thousands of regions */
    FILE **files = (FILE **)calloc(nthreads, sizeof(FILE *));
    for (int i = 0; i < nthreads; ++i) {
        files[i] = tmpfile();
        if (!files[i]) {
            fprintf(stderr, "[countmut] error: cannot create temp output file\n");
            for (int j = 0; j < i; ++j) fclose(files[j]);
            free(files); bam_hdr_destroy(hdr);
            if (fp != stdout) fclose(fp);
            return 1;
        }
    }
    worker_t *workers = (worker_t *)calloc(nthreads, sizeof(worker_t));
    for (int i = 0; i < nthreads; ++i)
        worker_init(&workers[i], bam, fa, cfg->pad, cfg->bedfile, cfg->exclude,
                    cfg->read_expr, cfg->pile_expr, cfg->output_expr,
                    use_pool ? &tpool : NULL);
    for (int i = 0; i < nthreads; ++i) {
        if (!workers[i].fp || !workers[i].idx) {
            fprintf(stderr, "[countmut] error: cannot open BAM/index '%s'\n", bam);
            goto fail_workers;
        }
        if (!workers[i].fai) {
            fprintf(stderr, "[countmut] error: cannot load reference FASTA '%s'\n", fa);
            goto fail_workers;
        }
    }
    { /* success path continues below */ }
    goto input_ok;
fail_workers:
    for (int i = 0; i < nthreads; ++i) worker_free(&workers[i]);
    free(workers);
    for (int i = 0; i < nthreads; ++i) fclose(files[i]);
    free(files);
    if (use_pool && tpool.pool) hts_tpool_destroy(tpool.pool);
    bam_hdr_destroy(hdr);
    if (fp != stdout) fclose(fp);
    return 3;
input_ok:
    ;

    work_t s;
    s.cfg = cfg; s.hdr = hdr; s.bam = bam;
    s.regions = regs; s.nregions = nregions; s.done = 0; s.next = 0;
    s.spans = (span_t *)calloc(nregions, sizeof(span_t));
    s.files = files; s.workers = workers; s.nthreads = nthreads;

    pthread_t *tds = (pthread_t *)calloc(nthreads, sizeof(pthread_t));
    targ_t *targs = (targ_t *)calloc(nthreads, sizeof(targ_t));
    for (int i = 0; i < nthreads; ++i) {
        targs[i].s = &s; targs[i].wi = i;
        pthread_create(&tds[i], NULL, thread_main, &targs[i]);
    }
    for (int i = 0; i < nthreads; ++i) pthread_join(tds[i], NULL);
    free(targs);

    /* re-assemble worker temp files in GLOBAL region order using the recorded
     * spans (the dynamic queue claims regions out of order) */
    for (int i = 0; i < nthreads; ++i) fflush(files[i]);
    char buf[16384];
    for (int i = 0; i < nregions; ++i) {
        int w = s.spans[i].worker;
        long remain = s.spans[i].len;
        if (remain <= 0) continue;
        fseek(files[w], s.spans[i].off, SEEK_SET);
        while (remain > 0) {
            size_t chunk = remain > (long)sizeof(buf) ? sizeof(buf) : (size_t)remain;
            size_t got = fread(buf, 1, chunk, files[w]);
            if (got == 0) break;
            fwrite(buf, 1, got, fp);
            remain -= (long)got;
        }
    }
    for (int i = 0; i < nthreads; ++i) fclose(files[i]);

    free(s.spans); free(tds); free(files);
    for (int i = 0; i < nthreads; ++i) worker_free(&workers[i]);
    free(workers); free(regs);
    if (use_pool && tpool.pool) hts_tpool_destroy(tpool.pool);
    bam_hdr_destroy(hdr);
    if (fp != stdout) fclose(fp);
    if (is_sam) {   /* clean up the transcoded temp BAM + its index */
        char bai[1200];
        snprintf(bai, sizeof(bai), "%s.bai", sam_tmp);
        unlink(sam_tmp);
        unlink(bai);
    }
    return 0;
}
