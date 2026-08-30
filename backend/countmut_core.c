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

#include "sam.h"
#include "faidx.h"
#include "ksort.h"
#include "khash.h"
#include "countmut_core.h"

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

static char comp_base(char c) {
    switch (toupper((unsigned char)c)) {
    case 'A': return 'T';
    case 'T':
    case 'U': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    default: return 'N';
    }
}

static char *revcomp(const char *s, int len) {
    char *out = (char *)malloc(len + 1);
    for (int i = 0; i < len; ++i) out[i] = comp_base(s[len - 1 - i]);
    out[len] = 0;
    return out;
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
    if (cfg->req_flags && (b->core.flag & cfg->req_flags) != cfg->req_flags) return 1;
    if (cfg->excl_flags && (b->core.flag & cfg->excl_flags)) return 1;
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
} worker_t;

static void worker_init(worker_t *w, const char *bam, const char *fa, int pad,
                        const char *bedfile, const char *exclude) {
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
}

typedef struct { int tid, beg, end; } region_t;

typedef struct {
    const cm_config *cfg;
    bam_hdr_t *hdr;
    const char *bam;
    region_t *regions;
    int nregions;
    volatile int next;
    FILE **files;
    worker_t *workers;
    int nthreads;
} work_t;

/* Write the header line for the configured mode. */
static void write_header(FILE *fp, const cm_config *cfg) {
    if (cfg->mode == CM_MODE_MUTATION) {
        fputs("chrom\tpos\tstrand\tmotif\tu0\tu1\tu2\tm0\tm1\tm2", fp);
        if (cfg->save_rest) fputs("\to0\to1\to2", fp);
        fputc('\n', fp);
    } else if (cfg->mode == CM_MODE_BASE) {
        if (cfg->split_strand) fputs("chrom\tpos\tstrand\tref\tdepth\ta\tc\tg\tt\tn", fp);
        else fputs("chrom\tpos\tref\tdepth\ta\tc\tg\tt\tn", fp);
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

    int pos;
    site_t site;
    while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) {
        if (pos < beg || pos >= end) continue;
        if (w->last_tid != tid) {
            if (w->chr_seq) free(w->chr_seq);
            w->chr_seq = fai_fetch(w->fai, hdr->target_name[tid], &w->chr_len);
            w->last_tid = tid;
        }
        if (w->chr_len == 0 || pos >= w->chr_len) continue;
        /* BED region restriction (pbr -b include / -x exclude) */
        if (w->inc_bed && !bed_overlap(w->inc_bed, hdr->target_name[tid], pos, pos + 1)) continue;
        if (w->exc_bed && bed_overlap(w->exc_bed, hdr->target_name[tid], pos, pos + 1)) continue;
        const int n = n_plp[0];
        if (n == 0) continue;
        const char ref_ch = (char)toupper((unsigned char)w->chr_seq[pos]);

        if (cfg->mode == CM_MODE_MUTATION && cfg->ref_base && toupper((unsigned char)ref_ch) != cfg->ref_base)
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
            if (read_fails(cfg, b)) { site.fail[s]++; continue; }
            if (p->is_refskip) { site.refskip[s]++; continue; }
            if (p->is_del) { site.del[s]++; continue; }
            if (p->qpos < 0 || p->qpos >= b->core.l_qseq) continue;
            int qpos = p->qpos, qlen = b->core.l_qseq;
            if (s == 0) { if (qpos < cfg->trim_start || qlen - qpos <= cfg->trim_end) continue; }
            else { if (qpos < cfg->trim_end || qlen - qpos <= cfg->trim_start) continue; }
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
        if (cfg->mode == CM_MODE_MUTATION) {
            int refi = base_to_index((char)cfg->ref_base), muti = base_to_index((char)cfg->mut_base);
            int mlen = cfg->pad * 2 + 1;
            for (int s = 0; s < 2; ++s) {
                if (s == 0 && !emit_plus) continue;
                if (s == 1 && !emit_minus) continue;
                int u0 = site.cnt[s][0][refi], m0 = site.cnt[s][0][muti];
                int u1 = site.cnt[s][1][refi], m1 = site.cnt[s][1][muti];
                int u2 = site.cnt[s][2][refi], m2 = site.cnt[s][2][muti];
                int o0 = site.cnt[s][0][0]+site.cnt[s][0][1]+site.cnt[s][0][2]+site.cnt[s][0][3]+site.cnt[s][0][4]-u0-m0;
                int o1 = site.cnt[s][1][0]+site.cnt[s][1][1]+site.cnt[s][1][2]+site.cnt[s][1][3]+site.cnt[s][1][4]-u1-m1;
                int o2 = site.cnt[s][2][0]+site.cnt[s][2][1]+site.cnt[s][2][2]+site.cnt[s][2][3]+site.cnt[s][2][4]-u2-m2;
                if (u1 + m1 + u2 + m2 <= 0) continue;
                for (int k2 = pos - cfg->pad; k2 < pos + cfg->pad + 1; ++k2) {
                    w->motif_buf[k2 - (pos - cfg->pad)] =
                        (k2 < 0 || k2 >= w->chr_len) ? 'N' : (char)toupper((unsigned char)w->chr_seq[k2]);
                }
                w->motif_buf[mlen] = 0;
                if (s == 1) {
                    char *rc = revcomp(w->motif_buf, mlen);
                    fprintf(fp, "%s\t%d\t-\t%s\t%d\t%d\t%d\t%d\t%d\t%d",
                            hdr->target_name[tid], pos + 1, rc, u0, u1, u2, m0, m1, m2);
                    free(rc);
                } else {
                    fprintf(fp, "%s\t%d\t+\t%s\t%d\t%d\t%d\t%d\t%d\t%d",
                            hdr->target_name[tid], pos + 1, w->motif_buf, u0, u1, u2, m0, m1, m2);
                }
                if (cfg->save_rest) fprintf(fp, "\t%d\t%d\t%d", o0, o1, o2);
                fputc('\n', fp);
            }
        } else if (cfg->mode == CM_MODE_BASE) {
            if (cfg->split_strand) {
                for (int s = 0; s < 2; ++s) {
                    if (s == 0 && !emit_plus) continue;
                    if (s == 1 && !emit_minus) continue;
                    int dep = site.cnt[s][0][0]+site.cnt[s][0][1]+site.cnt[s][0][2]+site.cnt[s][0][3]+site.cnt[s][0][4];
                    int t_ins = site.ins[s], t_del = site.del[s], t_rs = site.refskip[s], t_fl = site.fail[s];
                    if (dep + t_rs + t_del + t_ins + t_fl == 0) continue;
                    fprintf(fp, "%s\t%d\t%c\t%c\t%d\t%d\t%d\t%d\t%d\t%d",
                            hdr->target_name[tid], pos + 1, s ? '-' : '+', ref_ch, dep,
                            site.cnt[s][0][0], site.cnt[s][0][1], site.cnt[s][0][2],
                            site.cnt[s][0][3], site.cnt[s][0][4]);
                    if (cfg->count_indels) fprintf(fp, "\t%d\t%d\t%d\t%d", t_ins, t_del, t_rs, t_fl);
                    fputc('\n', fp);
                }
            } else {
                int cnt[5] = {0}; int t_ins = 0, t_del = 0, t_rs = 0, t_fl = 0;
                if (emit_plus) {
                    for (int b = 0; b < 5; ++b) cnt[b] += site.cnt[0][0][b];
                    t_ins += site.ins[0]; t_del += site.del[0]; t_rs += site.refskip[0]; t_fl += site.fail[0];
                }
                if (emit_minus) {
                    for (int b = 0; b < 5; ++b) cnt[b] += site.cnt[1][0][b];
                    t_ins += site.ins[1]; t_del += site.del[1]; t_rs += site.refskip[1]; t_fl += site.fail[1];
                }
                int dep = cnt[0]+cnt[1]+cnt[2]+cnt[3]+cnt[4];
                if (dep + t_rs + t_del + t_ins + t_fl == 0) continue;
                if (cfg->min_depth > 0 && dep < cfg->min_depth) continue;
                fprintf(fp, "%s\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d",
                        hdr->target_name[tid], pos + 1, ref_ch, dep, cnt[0], cnt[1], cnt[2], cnt[3], cnt[4]);
                if (cfg->count_indels) fprintf(fp, "\t%d\t%d\t%d\t%d", t_ins, t_del, t_rs, t_fl);
                fputc('\n', fp);
            }
        } else { /* allele */
            int cnt[5] = {0};
            if (emit_plus) for (int b = 0; b < 5; ++b) cnt[b] += site.cnt[0][0][b];
            if (emit_minus) for (int b = 0; b < 5; ++b) cnt[b] += site.cnt[1][0][b];
            int dep = cnt[0]+cnt[1]+cnt[2]+cnt[3]+cnt[4];
            if (dep <= 0) continue;
            if (cfg->min_depth > 0 && dep < cfg->min_depth) continue;
            int refi = base_to_index(ref_ch), ref_n = cnt[refi], best = -1, bn = 0;
            for (int i = 0; i < 4; ++i) if (i != refi && cnt[i] > bn) { bn = cnt[i]; best = i; }
            if (cfg->vcf) {
                if (best < 0 || bn < cfg->min_allele_support) continue;
                const char *alts = "ACGT";
                fprintf(fp, "%s\t%d\t.\t%c\t%c\t.\tPASS\t.\tGT:AD\t0/1:%d,%d\n",
                        hdr->target_name[tid], pos + 1, ref_ch, alts[best], ref_n, bn);
            } else {
                fprintf(fp, "%s\t%d\t%c\t%d\t%d\t%c\t%d\n",
                        hdr->target_name[tid], pos + 1, ref_ch, dep, ref_n,
                        best < 0 ? '.' : "ACGT"[best], best < 0 ? 0 : bn);
            }
        }
    }
    free(n_plp); free(plp);
    bam_mplp_destroy(mplp);
    if (aux.itr) bam_itr_destroy(aux.itr);
}

typedef struct { work_t *s; int wi; } targ_t;

static void *thread_main(void *arg) {
    targ_t *ta = (targ_t *)arg;
    worker_t *w = &ta->s->workers[ta->wi];
    work_t *s = ta->s;
    for (;;) {
        int i = __sync_fetch_and_add(&s->next, 1);
        if (i >= s->nregions) break;
        count_interval(w, s->cfg, s->hdr, s->files[i], s->regions[i].tid,
                       s->regions[i].beg, s->regions[i].end);
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

int cm_run(const cm_config *cfg, const char *bam, const char *fa, const char *out_path, const char *region) {
    FILE *fp = (out_path && strcmp(out_path, "-") != 0) ? fopen(out_path, "w") : stdout;
    if (!fp) return 1;
    BGZF *hfp = bgzf_open(bam, "r");
    bam_hdr_t *hdr = bam_hdr_read(hfp);
    bgzf_close(hfp);
    if (!hdr) { if (fp != stdout) fclose(fp); return 3; }

    int nthreads = cfg->threads < 1 ? 1 : cfg->threads;
    int nregions = 0;
    region_t *regs = build_regions(hdr, nthreads, region, &nregions);
    if (!regs) { bam_hdr_destroy(hdr); if (fp != stdout) fclose(fp); return 4; }
    if (nregions > 1 && nthreads > nregions) nthreads = nregions;

    write_header(fp, cfg);

    FILE **files = (FILE **)calloc(nregions, sizeof(FILE *));
    for (int i = 0; i < nregions; ++i) files[i] = tmpfile();
    worker_t *workers = (worker_t *)calloc(nthreads, sizeof(worker_t));
    for (int i = 0; i < nthreads; ++i)
        worker_init(&workers[i], bam, fa, cfg->pad, cfg->bedfile, cfg->exclude);

    work_t s;
    s.cfg = cfg; s.hdr = hdr; s.bam = bam;
    s.regions = regs; s.nregions = nregions; s.next = 0;
    s.files = files; s.workers = workers; s.nthreads = nthreads;

    pthread_t *tds = (pthread_t *)calloc(nthreads, sizeof(pthread_t));
    targ_t *targs = (targ_t *)calloc(nthreads, sizeof(targ_t));
    for (int i = 0; i < nthreads; ++i) {
        targs[i].s = &s; targs[i].wi = i;
        pthread_create(&tds[i], NULL, thread_main, &targs[i]);
    }
    for (int i = 0; i < nthreads; ++i) pthread_join(tds[i], NULL);
    free(targs);

    for (int i = 0; i < nregions; ++i) {
        fflush(files[i]); rewind(files[i]);
        char buf[16384]; size_t r;
        while ((r = fread(buf, 1, sizeof(buf), files[i])) > 0) fwrite(buf, 1, r, fp);
        fclose(files[i]);
    }

    free(tds); free(files);
    for (int i = 0; i < nthreads; ++i) worker_free(&workers[i]);
    free(workers); free(regs);
    bam_hdr_destroy(hdr);
    if (fp != stdout) fclose(fp);
    return 0;
}
