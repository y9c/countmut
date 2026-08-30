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

#include "sam.h"
#include "faidx.h"
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
static int read_trim_skip(const bam1_t *b, int qpos, int r1_end, int r2_start) {
    if (!(b->core.flag & BAM_FPAIRED)) return 0;
    if (r1_end > 0 && (b->core.flag & BAM_FREAD1))
        return qpos >= (int)b->core.l_qseq - r1_end;
    if (r2_start > 0 && (b->core.flag & BAM_FREAD2))
        return qpos < r2_start;
    return 0;
}

/* per-site accumulator */
typedef struct {
    int cnt[2][3][5]; /* [strand][category][base] */
    int ins[2], del[2], refskip[2], fail[2];
} site_t;

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

/* read-level filters (samtools reqflags/exclflags + mapq + NS + bisulfite tags) */
static int read_fails(const cm_config *cfg, const bam1_t *b) {
    if (cfg->req_flags && (b->core.flag & (uint32_t)cfg->req_flags) != (uint32_t)cfg->req_flags) return 1;
    if (cfg->excl_flags && (b->core.flag & (uint32_t)cfg->excl_flags)) return 1;
    if ((int)b->core.qual < cfg->min_mapq) return 1;
    if (cfg->max_sub >= 0) {
        uint8_t *aux = bam_aux_get(b, "NS");
        if (aux && bam_aux2i(aux) > cfg->max_sub) return 1;
    }
    if (cfg->max_unc >= 0 && cfg->min_con >= 0) {
        uint8_t *zf = bam_aux_get(b, "Zf");
        uint8_t *yf = bam_aux_get(b, "Yf");
        if (zf && yf) {
            if (!(bam_aux2i(zf) <= cfg->max_unc && bam_aux2i(yf) >= cfg->min_con)) return 1;
        }
    }
    return 0;
}

typedef struct {
    BGZF *fp;
    hts_itr_t *itr;
    int beg, end;
} aux_t;

static int read_bam(void *data, bam1_t *b) {
    aux_t *aux = (aux_t *)data;
    int ret = aux->itr ? bam_itr_next(aux->fp, aux->itr, b) : bam_read1(aux->fp, b);
    return ret;
}

KHASH_INIT(qn, char *, int, 1, kh_str_hash_func, kh_str_hash_equal)

/* read-filter memo: (ref_start, qname) -> pass/fail, so read_fails() is
 * computed once per read instead of once per pileup position. */
typedef struct { int64_t pos; const char *qn; } rf_key;
static inline khint_t rf_hash(rf_key k) {
    return kh_int64_hash_func(k.pos) ^ kh_str_hash_func(k.qn);
}
static inline int rf_equal(rf_key a, rf_key b) {
    return a.pos == b.pos && strcmp(a.qn, b.qn) == 0;
}
KHASH_INIT(rfc, rf_key, int, 1, rf_hash, rf_equal)
#define RF_CAP (1 << 15)

/* -e read-constant memo for the pileup engine, keyed by the pileup slot
 * pointer (p->b): htslib keeps the same bam1_t for one read across all the
 * positions it covers, so an int-keyed slot hash lets us evaluate a
 * read-constant expression once per read and reuse it for every appearance.
 * A pos/qlen/qname verify guards against a recycled buffer now holding a
 * different read.  Unlike the RF_CAP-based cache this survives deep hotspots
 * (no global eviction: each mplp slot holds at most one read at a time). */
typedef struct { int64_t pos; int qlen; char *qn; int pass; } expr_cc_t;
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
    int *sel, *mapq_a, *r1_a, *q_a, sel_cap;
    char *motif_buf;
    char *chr_seq; int chr_len, last_tid;
    BGZF *fp; hts_idx_t *idx; faidx_t *fai;
    void *inc_bed, *exc_bed;
    cm_expr *expr;   /* Lua -e / -p filters (NULL when none) */
    khash_t(rfc) *rfc; /* read_fails memo (pileup engine) */
    khash_t(pex) *pexc; /* -e read-constant memo keyed by pileup slot (pileup) */
} worker_t;

static void worker_init(worker_t *w, const char *bam, const char *fa, int pad,
                        const char *bedfile, const char *exclude,
                        const char *read_expr, const char *pile_expr) {
    w->kh = kh_init(qn);
    w->sel = w->mapq_a = w->r1_a = w->q_a = NULL;
    w->sel_cap = 0;
    w->motif_buf = (char *)malloc(2 * pad + 2);
    w->chr_seq = NULL; w->chr_len = 0; w->last_tid = -1;
    w->fp = bgzf_open(bam, "r");
    w->idx = bam_index_load(bam);
    w->fai = fai_load(fa);
    w->inc_bed = bedfile ? bed_read(bedfile) : NULL;
    w->exc_bed = exclude ? bed_read(exclude) : NULL;
    w->expr = cm_expr_new(read_expr, pile_expr);
    w->rfc = kh_init(rfc);
    w->pexc = kh_init(pex);
}

static void worker_free(worker_t *w) {
    if (w->sel) free(w->sel);
    if (w->mapq_a) free(w->mapq_a);
    if (w->r1_a) free(w->r1_a);
    if (w->q_a) free(w->q_a);
    if (w->motif_buf) free(w->motif_buf);
    if (w->chr_seq) free(w->chr_seq);
    kh_destroy(qn, w->kh);
    if (w->exc_bed) bed_destroy(w->exc_bed);
    if (w->inc_bed) bed_destroy(w->inc_bed);
    if (w->fai) fai_destroy(w->fai);
    if (w->idx) hts_idx_destroy(w->idx);
    if (w->fp) bgzf_close(w->fp);
    cm_expr_free(w->expr);
    if (w->rfc) {
        for (khint_t k = kh_begin(w->rfc); k != kh_end(w->rfc); ++k)
            if (kh_exist(w->rfc, k)) free((void *)(uintptr_t)kh_key(w->rfc, k).qn);
        kh_destroy(rfc, w->rfc);
    }
    if (w->pexc) {
        for (khint_t k = kh_begin(w->pexc); k != kh_end(w->pexc); ++k)
            if (kh_exist(w->pexc, k)) free(kh_val(w->pexc, k).qn);
        kh_destroy(pex, w->pexc);
    }
}

/* read_fails() memoized by (ref_start, qname).  Only used when an aux-tag
 * filter (NS / Zf / Yf) is active -- otherwise read_fails() is cheap and the
 * cache would only add overhead.  The cached value is exactly read_fails() for
 * that same read, so input->output semantics are unchanged. */
#define read_fails_cached(w, cfg, b) \
    (((cfg)->max_sub < 0 && (cfg)->max_unc < 0 && (cfg)->min_con < 0) ? \
        read_fails(cfg, b) : _read_fails_cached(w, cfg, b))

static int _read_fails_cached(worker_t *w, const cm_config *cfg, const bam1_t *b) {
    rf_key key = { b->core.pos, bam_get_qname(b) };
    khint_t k = kh_get(rfc, w->rfc, key);
    if (k != kh_end(w->rfc)) return kh_val(w->rfc, k);
    int fails = read_fails(cfg, b);
    if (kh_size(w->rfc) >= RF_CAP) {
        for (khint_t it = kh_begin(w->rfc); it != kh_end(w->rfc); ++it)
            if (kh_exist(w->rfc, it)) free((void *)(uintptr_t)kh_key(w->rfc, it).qn);
        kh_clear(rfc, w->rfc);
    }
    int ret; k = kh_put(rfc, w->rfc, key, &ret);
    if (ret) kh_key(w->rfc, k).qn = strdup(bam_get_qname(b));
    kh_val(w->rfc, k) = fails;
    return fails;
}

/* -e read filter, memoized by (ref_start, qname) WHEN the expression is
 * read-constant (no per-base qpos/bq/base/ref/dist).  A read-constant
 * expression yields the same result at every base of a read, so caching makes
 * the pileup engine evaluate it once per read instead of once per position
 * (the read-walk engine already does once per read).  Per-base expressions are
 * evaluated at every aligned base, uncached. */
static int expr_pass(worker_t *w, const bam1_t *b, const char *rname,
                     const char *mrname, int s, int qpos, char ref_ch) {
    cm_expr *x = w->expr;
    if (x == NULL || !cm_expr_has_read(x)) return 1;
    if (!cm_expr_read_constant(x))
        return cm_expr_read(x, b, rname, mrname, qpos, s ? -1 : 1, ref_ch);
    /* read-constant: memoize by the pileup slot pointer (stable per read
     * across its span) with a pos/qlen/qname verify against recycling. */
    uintptr_t slot = (uintptr_t)(const void *)b;
    khint_t k = kh_get(pex, w->pexc, slot);
    if (k != kh_end(w->pexc)) {
        expr_cc_t *cc = &kh_val(w->pexc, k);
        if (cc->pos == b->core.pos && cc->qlen == (int)b->core.l_qseq
            && cc->qn && strcmp(cc->qn, bam_get_qname(b)) == 0)
            return cc->pass;
        free(cc->qn); cc->qn = NULL;
    }
    int pass = cm_expr_read(x, b, rname, mrname, 0, s ? -1 : 1, 'N');
    if (k == kh_end(w->pexc)) {
        int ret; k = kh_put(pex, w->pexc, slot, &ret);
        memset(&kh_val(w->pexc, k), 0, sizeof(expr_cc_t)); /* fresh slots are uninitialized */
    }
    expr_cc_t *cc = &kh_val(w->pexc, k);
    if (cc->qn == NULL) cc->qn = strdup(bam_get_qname(b));
    cc->pos = b->core.pos;
    cc->qlen = (int)b->core.l_qseq;
    cc->pass = pass;
    return pass;
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

/* Write the header line for the configured mode. */
static void write_header(FILE *fp, const cm_config *cfg) {
    if (cfg->mode == CM_MODE_MUTATION) {
        fputs("chrom\tpos\tstrand\tmotif\tu0\tu1\tu2\tm0\tm1\tm2", fp);
        if (cfg->save_rest) fputs("\to0\to1\to2", fp);
        fputs("\tmutation_rate", fp);
        fputc('\n', fp);
    } else if (cfg->mode == CM_MODE_BASE) {
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
    if (cfg->mode == CM_MODE_MUTATION) {
        int refi = base_to_index((char)cfg->ref_base), muti = base_to_index((char)cfg->mut_base);
        int mlen = cfg->pad * 2 + 1;
        int motif_ready = 0;   /* build the reference-forward window once */
        for (int s = 0; s < 2; ++s) {
            if (s == 0 && !emit_plus) continue;
            if (s == 1 && !emit_minus) continue;
            int u0 = site->cnt[s][0][refi], m0 = site->cnt[s][0][muti];
            int u1 = site->cnt[s][1][refi], m1 = site->cnt[s][1][muti];
            int u2 = site->cnt[s][2][refi], m2 = site->cnt[s][2][muti];
            int ut = u0 + u1 + u2, mt = m0 + m1 + m2;   /* total unconverted / converted */
            int o0 = site->cnt[s][0][0]+site->cnt[s][0][1]+site->cnt[s][0][2]+site->cnt[s][0][3]+site->cnt[s][0][4]-u0-m0;
            int o1 = site->cnt[s][1][0]+site->cnt[s][1][1]+site->cnt[s][1][2]+site->cnt[s][1][3]+site->cnt[s][1][4]-u1-m1;
            int o2 = site->cnt[s][2][0]+site->cnt[s][2][1]+site->cnt[s][2][2]+site->cnt[s][2][3]+site->cnt[s][2][4]-u2-m2;
            if (u1 + m1 + u2 + m2 <= 0) continue;
            if (!motif_ready) {
                /* Reference-forward window for BOTH strands (bases are counted
                 * reference-forward), built at most once per site. */
                for (int k2 = (int)pos - cfg->pad; k2 < (int)pos + cfg->pad + 1; ++k2) {
                    w->motif_buf[k2 - ((int)pos - cfg->pad)] =
                        (k2 < 0 || k2 >= w->chr_len) ? 'N'
                        : (char)toupper((unsigned char)w->chr_seq[k2]);
                }
                w->motif_buf[mlen] = 0;
                motif_ready = 1;
            }
            fprintf(fp, "%s\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d",
                    hdr->target_name[tid], (int)pos + 1, s ? '-' : '+', w->motif_buf,
                    u0, u1, u2, m0, m1, m2);
            if (cfg->save_rest) fprintf(fp, "\t%d\t%d\t%d", o0, o1, o2);
            fprintf(fp, "\t%.4f",
                    (ut + mt) ? (double)mt / (double)(ut + mt)
                              : strtod("nan", (char **)NULL));
            fputc('\n', fp);
        }
    } else if (cfg->mode == CM_MODE_BASE) {
        if (!cfg->strandless) {
            for (int s = 0; s < 2; ++s) {
                if (s == 0 && !emit_plus) continue;
                if (s == 1 && !emit_minus) continue;
                int dep = site->cnt[s][0][0]+site->cnt[s][0][1]+site->cnt[s][0][2]+site->cnt[s][0][3]+site->cnt[s][0][4];
                int t_ins = site->ins[s], t_del = site->del[s], t_rs = site->refskip[s], t_fl = site->fail[s];
                if (dep + t_rs + t_del + t_ins + t_fl == 0) continue;
                if (cfg->min_depth > 0 && dep < cfg->min_depth) continue;
                fprintf(fp, "%s\t%d\t%c\t%c\t%d\t%d\t%d\t%d\t%d\t%d",
                        hdr->target_name[tid], (int)pos + 1, s ? '-' : '+', ref_ch, dep,
                        site->cnt[s][0][0], site->cnt[s][0][1], site->cnt[s][0][2],
                        site->cnt[s][0][3], site->cnt[s][0][4]);
                if (cfg->count_indels) fprintf(fp, "\t%d\t%d\t%d\t%d", t_ins, t_del, t_rs, t_fl);
                fputc('\n', fp);
            }
        } else {
            int cnt[5] = {0}; int t_ins = 0, t_del = 0, t_rs = 0, t_fl = 0;
            if (emit_plus) {
                for (int b = 0; b < 5; ++b) cnt[b] += site->cnt[0][0][b];
                t_ins += site->ins[0]; t_del += site->del[0]; t_rs += site->refskip[0]; t_fl += site->fail[0];
            }
            if (emit_minus) {
                for (int b = 0; b < 5; ++b) cnt[b] += site->cnt[1][0][b];
                t_ins += site->ins[1]; t_del += site->del[1]; t_rs += site->refskip[1]; t_fl += site->fail[1];
            }
            int dep = cnt[0]+cnt[1]+cnt[2]+cnt[3]+cnt[4];
            if (dep + t_rs + t_del + t_ins + t_fl == 0) return;
            if (cfg->min_depth > 0 && dep < cfg->min_depth) return;
            fprintf(fp, "%s\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d",
                    hdr->target_name[tid], (int)pos + 1, ref_ch, dep, cnt[0], cnt[1], cnt[2], cnt[3], cnt[4]);
            if (cfg->count_indels) fprintf(fp, "\t%d\t%d\t%d\t%d", t_ins, t_del, t_rs, t_fl);
            fputc('\n', fp);
        }
    } else { /* allele */
        int cnt[5] = {0};
        if (emit_plus) for (int b = 0; b < 5; ++b) cnt[b] += site->cnt[0][0][b];
        if (emit_minus) for (int b = 0; b < 5; ++b) cnt[b] += site->cnt[1][0][b];
        int dep = cnt[0]+cnt[1]+cnt[2]+cnt[3]+cnt[4];
        if (dep <= 0) return;
        if (cfg->min_depth > 0 && dep < cfg->min_depth) return;
        int refi = base_to_index(ref_ch), ref_n = cnt[refi], best = -1, bn = 0;
        for (int i = 0; i < 4; ++i) if (i != refi && cnt[i] > bn) { bn = cnt[i]; best = i; }
        if (cfg->vcf) {
            if (best < 0 || bn < cfg->min_allele_support) return;
            const char *alts = "ACGT";
            fprintf(fp, "%s\t%d\t.\t%c\t%c\t.\tPASS\t.\tGT:AD\t0/1:%d,%d\n",
                    hdr->target_name[tid], (int)pos + 1, ref_ch, alts[best], ref_n, bn);
        } else {
            if (bn < cfg->min_allele_support) { best = -1; bn = 0; }
            fprintf(fp, "%s\t%d\t%c\t%d\t%d\t%c\t%d\n",
                    hdr->target_name[tid], (int)pos + 1, ref_ch, dep, ref_n,
                    best < 0 ? '.' : "ACGT"[best], best < 0 ? 0 : bn);
        }
    }
}

/* Evaluate the -p (pile/site) filter for a fully-built site_t.  Computes the
 * A/C/G/T/N totals (both strands, all quality tiers) plus indels, and the
 * reference window for mutation mode.  Returns 1 = keep, 0 = omit. */
static int expr_pile_apply(cm_expr *x, const cm_config *cfg, worker_t *w,
                           const site_t *site, int64_t pos, char ref_ch) {
    int cnt[5] = {0};
    int ins = 0, del = 0, rs = 0, fl = 0;
    for (int s = 0; s < 2; ++s) {
        for (int c = 0; c < 3; ++c)
            for (int b = 0; b < 5; ++b) cnt[b] += site->cnt[s][c][b];
        ins += site->ins[s]; del += site->del[s]; rs += site->refskip[s]; fl += site->fail[s];
    }
    const char *motif = NULL;
    if (cfg->mode == CM_MODE_MUTATION && w->chr_len > 0) {
        int mlen = cfg->pad * 2 + 1;
        for (int k2 = (int)pos - cfg->pad; k2 < (int)pos + cfg->pad + 1; ++k2) {
            w->motif_buf[k2 - ((int)pos - cfg->pad)] =
                (k2 < 0 || k2 >= w->chr_len) ? 'N' : (char)toupper((unsigned char)w->chr_seq[k2]);
        }
        w->motif_buf[mlen] = 0;
        motif = w->motif_buf;
    }
    return cm_expr_pile(x, pos, ref_ch, motif, cnt, ins, del, rs, fl);
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
    aux.fp = w->fp;
    aux.beg = beg; aux.end = end;
    aux.itr = w->idx ? bam_itr_queryi(w->idx, tid, beg, end) : NULL;

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

        if (cfg->mode == CM_MODE_MUTATION && cfg->ref_base && ref_ch != cfg->ref_base)
            continue;

        if (n > w->sel_cap) {
            w->sel = (int *)realloc(w->sel, n * sizeof(int));
            w->mapq_a = (int *)realloc(w->mapq_a, n * sizeof(int));
            w->r1_a = (int *)realloc(w->r1_a, n * sizeof(int));
            w->q_a = (int *)realloc(w->q_a, n * sizeof(int));
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
            if (read_fails_cached(w, cfg, b)) { site.fail[s]++; continue; }
            if (p->is_refskip) { site.refskip[s]++; continue; }
            if (p->is_del) { site.del[s]++; continue; }
            if (p->qpos < 0 || p->qpos >= b->core.l_qseq) continue;
            int qpos = p->qpos, qlen = b->core.l_qseq;
            if (s == 0) { if (qpos < cfg->trim_fragment_start || qlen - qpos <= cfg->trim_fragment_end) continue; }
            else { if (qpos < cfg->trim_fragment_end || qlen - qpos <= cfg->trim_fragment_start) continue; }
            if (read_trim_skip(b, qpos, cfg->trim_r1_end, cfg->trim_r2_start)) continue;
            /* -e read filter (once per read when read-constant via exprc memo,
             * else per aligned base; same spot as the Python engine) */
            if (!expr_pass(w, b, hdr->target_name[tid],
                           (b->core.mtid >= 0 && b->core.mtid < hdr->n_targets)
                               ? hdr->target_name[b->core.mtid] : "",
                           s, qpos, ref_ch))
                continue;
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
            int base_i = nt16_index(nt);   /* stored SEQ is reference-forward */
            int qual = (int)bam_get_qual(b)[p->qpos];
            if (cfg->mode == CM_MODE_MUTATION) {
                site.cnt[s][(qual >= cfg->min_baseq) ? 2 : 0][base_i]++;
            } else {
                site.cnt[s][0][base_i]++;
            }
        }

        /* ---------- emit ---------- */
        const int emit_plus = cfg->strand_process != CM_STRAND_REVERSE;
        const int emit_minus = cfg->strand_process != CM_STRAND_FORWARD;
        /* -p site filter */
        if (w->expr && cm_expr_has_pile(w->expr)
            && !expr_pile_apply(w->expr, cfg, w, &site, pos, ref_ch))
            continue;
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
static inline khint_t posq_hash(posq_key k) {
    return kh_int64_hash_func(k.pos) ^ (khint_t)k.qid;
}
static inline int posq_equal(posq_key a, posq_key b) {
    return a.pos == b.pos && a.qid == b.qid;
}
KHASH_INIT(posq, posq_key, int, 1, posq_hash, posq_equal)
KHASH_INIT(qn2id, char *, int, 1, kh_str_hash_func, kh_str_hash_equal)

/* pos -> sitemap slot */
KHASH_INIT(posi, khint64_t, int, 1, kh_int64_hash_func, kh_int64_hash_equal)

/* winner of a (pos,qname) dedup bucket */
typedef struct { int mapq, r1, qual, strand, base; } rw_w;

/* growable map pos -> site_t (one entry per visited reference position) */
typedef struct {
    khash_t(posi) *pm;
    int64_t *spos;
    site_t *st;
    int n, cap;
} sitemap_t;

static void sm_init(sitemap_t *m) {
    m->pm = kh_init(posi); m->spos = NULL; m->st = NULL; m->n = m->cap = 0;
}
static void sm_free(sitemap_t *m) {
    kh_destroy(posi, m->pm); free(m->spos); free(m->st);
    memset(m, 0, sizeof(*m));
}
/* NOTE: the returned pointer is only valid until the next sm_get() (a later
 * insert may realloc), so it must be used immediately and never retained. */
static site_t *sm_get(sitemap_t *m, int64_t pos) {
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

static uint32_t cigar_ref_len(const bam1_t *b) {
    const uint32_t *cig = bam_get_cigar(b);
    uint32_t r = 0;
    for (int i = 0; i < b->core.n_cigar; ++i) {
        int op = (int)bam_cigar_op_p(&cig[i]);
        if (op == 0 || op == 2 || op == 3 || op == 7 || op == 8)
            r += (uint32_t)bam_cigar_oplen_p(&cig[i]);
    }
    return r;
}
static int cigar_has_indels(const bam1_t *b) {
    const uint32_t *cig = bam_get_cigar(b);
    for (int i = 0; i < b->core.n_cigar; ++i) {
        int op = (int)bam_cigar_op_p(&cig[i]);
        if (op == 1 || op == 2 || op == 3) return 1;
    }
    return 0;
}
static uint32_t cigar_leading_softclips(const bam1_t *b) {
    const uint32_t *cig = bam_get_cigar(b);
    uint32_t sc = 0;
    for (int i = 0; i < b->core.n_cigar; ++i) {
        int op = (int)bam_cigar_op_p(&cig[i]);
        if (op == 4) sc += (uint32_t)bam_cigar_oplen_p(&cig[i]);
        else break;   /* soft-clips only at the 5' end before the first M */
    }
    return sc;
}
static int tgt_lower_bound(const int *a, int n, int v) {
    int lo = 0, hi = n;
    while (lo < hi) { int mid = (lo + hi) >> 1; if (a[mid] < v) lo = mid + 1; else hi = mid; }
    return lo;
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
    if (!already_target && cfg->mode == CM_MODE_MUTATION && cfg->ref_base) {
        if (ref_pos >= w->chr_len || w->chr_seq[ref_pos] != cfg->ref_base) return wins;
    }
    if (s == 0) {
        if ((int)qpos < cfg->trim_fragment_start || (int)qlen - (int)qpos <= cfg->trim_fragment_end) return wins;
    } else {
        if ((int)qpos < cfg->trim_fragment_end || (int)qlen - (int)qpos <= cfg->trim_fragment_start) return wins;
    }
    if (read_trim_skip(b, (int)qpos, cfg->trim_r1_end, cfg->trim_r2_start)) return wins;
    if (w->expr && cm_expr_has_read(w->expr) && !cm_expr_read_constant(w->expr)
        && !cm_expr_read(w->expr, b, hdr->target_name[tid],
                         (b->core.mtid >= 0 && b->core.mtid < hdr->n_targets)
                             ? hdr->target_name[b->core.mtid] : "",
                         (int)qpos, s ? -1 : 1,
                         (ref_pos >= 0 && ref_pos < w->chr_len) ? w->chr_seq[ref_pos] : 'N'))
        return wins;
    uint8_t nt = bam_seqi(bam_get_seq(b), qpos);
    int base_i = nt16_index(nt);   /* stored SEQ is reference-forward */
    int qual = (int)bam_get_qual(b)[qpos];
    if (direct) {
        int cat = (cfg->mode == CM_MODE_MUTATION)
            ? ((qual >= cfg->min_baseq) ? 2 : 0) : 0;
        sm_get(sm, ref_pos)->cnt[s][cat][base_i]++;
        return wins;
    }
    int mapq = (int)b->core.qual;
    int r1 = (b->core.flag & BAM_FREAD1) ? 1 : 0;
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
        wins[idx].strand = s; wins[idx].base = base_i;
        kh_val(h, kh) = idx;
    } else {
        int j = kh_val(h, kh);
        if (better(mapq, r1, qual, wins[j].mapq, wins[j].r1, wins[j].qual)) {
            wins[j].mapq = mapq; wins[j].r1 = r1; wins[j].qual = qual;
            wins[j].strand = s; wins[j].base = base_i;
        }
    }
    return wins;
}

static void count_interval_readwalk(worker_t *w, const cm_config *cfg, bam_hdr_t *hdr,
                                    FILE *fp, int tid, int beg, int end) {
    aux_t aux;
    aux.fp = w->fp; aux.beg = beg; aux.end = end;
    aux.itr = w->idx ? bam_itr_queryi(w->idx, tid, beg, end) : NULL;

    if (w->last_tid != tid) load_chr_seq(w, hdr, tid);

    /* sorted target positions (mutation mode) for the indel-free fast path */
    int *tgt = NULL; int tgt_n = 0;
    if (cfg->mode == CM_MODE_MUTATION && cfg->ref_base && end > beg) {
        tgt = (int *)malloc((size_t)(end - beg) * sizeof(int));
        for (int p = beg; p < end && p < w->chr_len; ++p)
            if (w->chr_seq[p] == (char)cfg->ref_base)
                tgt[tgt_n++] = p;
    }

    khash_t(posq) *h = kh_init(posq);
    khash_t(qn2id) *qnids = kh_init(qn2id);
    int qname_n = 0;
    rw_w *wins = NULL; int wins_cap = 0, wins_n = 0;
    sitemap_t sm; sm_init(&sm);

    bam1_t *b = bam_init1();
    int ret;
    while ((ret = (aux.itr ? bam_itr_next(aux.fp, aux.itr, b) : bam_read1(aux.fp, b))) >= 0) {
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
                        if (cfg->mode == CM_MODE_MUTATION && cfg->ref_base) {
                            for (int64_t p = lo; p < hi; ++p)
                                if ((int)p < w->chr_len && w->chr_seq[p] == cfg->ref_base)
                                    sm_get(&sm, p)->fail[s]++;
                        } else {
                            for (int64_t p = lo; p < hi; ++p)
                                sm_get(&sm, p)->fail[s]++;
                        }
                    }
                    rcur += len;
                }
            }
            continue;
        }

        uint32_t qlen = b->core.l_qseq;
        if (qlen == 0) continue;
        const char *qnbuf = bam_get_qname(b);
        khint_t qk = kh_get(qn2id, qnids, qnbuf);
        int qid;
        if (qk == kh_end(qnids)) {
            int r; char *cp = strdup(qnbuf); qk = kh_put(qn2id, qnids, cp, &r);
            qid = qname_n++; kh_val(qnids, qk) = qid;
        } else {
            qid = kh_val(qnids, qk);
        }
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
            /* same-contig mate: hybrid, if the insert geometry is usable */
            int insv = (int)b->core.isize; if (insv < 0) insv = -insv;
            if (insv > 0 && insv < 2 * (int)qlen) {
                ovl = 1;
                if (b->core.flag & BAM_FREAD1) olo = insv - (int)qlen;
                else ohi = 2 * (int)qlen - insv;
            }
        }
        /* else (paired but mate position unknown: mpos<0/mtid<0) ovl stays -1
         * -> every base goes through the hash = exact (never direct) */
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

        if (tgt && !cigar_has_indels(b)) {
            /* fast path: indel-free read -> jump straight to target positions */
            int64_t r_start = b->core.pos;
            uint32_t sc = cigar_leading_softclips(b);
            int64_t r_end = r_start + (int64_t)cigar_ref_len(b);
            int64_t lo = r_start > beg ? r_start : beg;
            int64_t hi = r_end < end ? r_end : end;
            for (int ti = tgt_lower_bound(tgt, tgt_n, (int)lo); ti < tgt_n; ++ti) {
                int64_t ref_pos = tgt[ti];
                if (ref_pos >= hi) break;
                uint32_t qpos = (uint32_t)((ref_pos - r_start) + sc);
                wins = rw_add_base(w, cfg, hdr, tid, b, s, ref_pos, qpos,
                                   qid, 1, RW_DIRECT(qpos), &sm,
                                   h, wins, &wins_cap, &wins_n);
            }
        } else {
            uint32_t qcur = 0;
            int64_t rcur = b->core.pos;
            for (int i = 0; i < b->core.n_cigar; ++i) {
                int op = (int)bam_cigar_op_p(&cig[i]); int len = (int)bam_cigar_oplen_p(&cig[i]);
                switch (op) {
                case 0: case 7: case 8: /* M, =, X -- matched bases */
                    for (int k = 0; k < len; ++k) {
                        int64_t ref_pos = rcur + k;
                        if (ref_pos < beg || ref_pos >= end) continue;
                        uint32_t qpos = qcur + (uint32_t)k;
                        if (qpos >= qlen) break;
                        wins = rw_add_base(w, cfg, hdr, tid, b, s, ref_pos, qpos,
                                           qid, 0, RW_DIRECT(qpos), &sm,
                                           h, wins, &wins_cap, &wins_n);
                    }
                    qcur += (uint32_t)len; rcur += len;
                    break;
                case 1: case 4: qcur += (uint32_t)len; break; /* I, S consume query */
                case 2: case 3: { /* D, N -- deletion / ref-skip */
                    int is_del = (op == 2);
                    if (rcur < end && rcur + len > beg) {
                        int64_t lo2 = rcur > beg ? rcur : beg;
                        int64_t hi2 = rcur + len < end ? rcur + len : end;
                        if (cfg->mode == CM_MODE_MUTATION && cfg->ref_base) {
                            for (int64_t p = lo2; p < hi2; ++p)
                                if ((int)p < w->chr_len && w->chr_seq[p] == cfg->ref_base) {
                                    if (is_del) sm_get(&sm, p)->del[s]++; else sm_get(&sm, p)->refskip[s]++;
                                }
                        } else {
                            for (int64_t p = lo2; p < hi2; ++p)
                                if (is_del) sm_get(&sm, p)->del[s]++; else sm_get(&sm, p)->refskip[s]++;
                        }
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
        int cat = (cfg->mode == CM_MODE_MUTATION)
            ? ((win->qual >= cfg->min_baseq) ? 2 : 0) : 0;
        site_t *st = sm_get(&sm, pos);
        st->cnt[win->strand][cat][win->base]++;
    }

    /* sort sites by position and emit */
    site_ord *ord = (site_ord *)malloc((size_t)(sm.n ? sm.n : 1) * sizeof(*ord));
    for (int i = 0; i < sm.n; ++i) { ord[i].pos = sm.spos[i]; ord[i].idx = i; }
    qsort(ord, (size_t)sm.n, sizeof(*ord), cmp_site_ord);
    const int emit_plus = cfg->strand_process != CM_STRAND_REVERSE;
    const int emit_minus = cfg->strand_process != CM_STRAND_FORWARD;
    for (int i = 0; i < sm.n; ++i) {
        int64_t pos = ord[i].pos;
        if (pos < 0 || pos >= w->chr_len) continue;
        char ref_ch = w->chr_seq[pos];   /* pre-uppercased */
        if (cfg->mode == CM_MODE_MUTATION && cfg->ref_base && ref_ch != cfg->ref_base) continue;
        if (w->inc_bed && !bed_overlap(w->inc_bed, hdr->target_name[tid], (int)pos, (int)pos + 1)) continue;
        if (w->exc_bed && bed_overlap(w->exc_bed, hdr->target_name[tid], (int)pos, (int)pos + 1)) continue;
        /* -p site filter */
        if (w->expr && cm_expr_has_pile(w->expr)
            && !expr_pile_apply(w->expr, cfg, w, &sm.st[ord[i].idx], pos, ref_ch))
            continue;
        emit_site(w, cfg, hdr, fp, tid, (int)pos, ref_ch, &sm.st[ord[i].idx], emit_plus, emit_minus);
    }

    /* cleanup */
    for (khint_t k = kh_begin(qnids); k != kh_end(qnids); ++k)
        if (kh_exist(qnids, k)) free((void *)(uintptr_t)kh_key(qnids, k));
    kh_destroy(qn2id, qnids);
    kh_destroy(posq, h);
    free(wins);
    free(ord);
    free(tgt);
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
static int detect_input_format(const char *path, int *is_sam) {
    unsigned char magic[4] = {0};
    gzFile gz = gzopen(path, "rb");
    if (gz == NULL) return -1;
    int n = (int)gzread(gz, magic, 4);
    gzclose(gz);
    if (n >= 4 && memcmp(magic, "BAM\1", 4) == 0) { *is_sam = 0; return 0; }
    if (n >= 4 && memcmp(magic, "CRAM", 4) == 0)   { *is_sam = 0; return 2; }
    *is_sam = 1; return 1;   /* SAM text (or empty -> header parse fails later with a clear error) */
}

static int transcode_sam_to_bam(const char *sam, char *tmp_bam, size_t cap) {
    htsFile *in = hts_open(sam, "r", NULL);
    if (in == NULL || in->is_bin) {
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

    BGZF *out = bgzf_open(tmp_bam, "w");
    if (out == NULL) {
        fprintf(stderr, "[countmut] error: cannot write temp BAM '%s'\n", tmp_bam);
        bam_hdr_destroy(hdr); hts_close(in);
        unlink(tmp_bam);
        return -1;
    }
    bam_hdr_write(out, hdr);
    bam1_t *b = bam_init1();
    int nrec = 0;
    while (sam_read1(in, hdr, b) >= 0) {
        bam_write1(out, b);
        ++nrec;
    }
    bgzf_close(out);
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(in);
    /* let the subset's own reader-driven indexer build the BAI (the hand-built
     * hts_idx_push path proved unreliable here) */
    if (bam_index_build(tmp_bam, 0) != 0) {
        fprintf(stderr, "[countmut] error: cannot index temp BAM '%s'\n", tmp_bam);
        unlink(tmp_bam);
        return -1;
    }
    fprintf(stderr, "[countmut] input is SAM: converted %d records -> %s\n", nrec, tmp_bam);
    return 0;
}

int cm_run(const cm_config *cfg, const char *bam, const char *fa, const char *out_path, const char *region) {
    FILE *fp = (out_path && strcmp(out_path, "-") != 0) ? fopen(out_path, "w") : stdout;
    if (!fp) return 1;
    int is_sam = 0;
    if (detect_input_format(bam, &is_sam) == 2) {
        fprintf(stderr,
                "[countmut] error: CRAM input '%s' is not supported by this "
                "self-contained core (it has no CRAM codec).  Convert it first:\n"
                "    samtools view -b %s -o out.bam\n"
                "(For CRAM with an embedded reference, samtools view also works "
                "without a separate FASTA.)\n", bam, bam);
        if (fp != stdout) fclose(fp);
        return 3;
    }
    char sam_tmp[1100] = {0};
    if (is_sam) {
        if (transcode_sam_to_bam(bam, sam_tmp, sizeof(sam_tmp)) != 0) {
            if (fp != stdout) fclose(fp);
            return 3;
        }
        bam = sam_tmp;   /* the rest of the run operates on the temp BAM */
    }
    BGZF *hfp = bgzf_open(bam, "r");
    if (!hfp) {
        fprintf(stderr, "[countmut] error: cannot open BAM file '%s'\n", bam);
        if (fp != stdout) fclose(fp);
        return 3;
    }
    bam_hdr_t *hdr = bam_hdr_read(hfp);
    bgzf_close(hfp);
    if (!hdr) {
        fprintf(stderr, "[countmut] error: cannot read BAM header from '%s'\n", bam);
        if (fp != stdout) fclose(fp);
        return 3;
    }

    int nthreads = cfg->threads < 1 ? 1 : cfg->threads;
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
                    cfg->read_expr, cfg->pile_expr);
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
