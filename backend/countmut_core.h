/* countmut_core.h -- C backend API for the unified countmut tool.
 *
 * The heavy lifting (reading BAM, pileup walk, base/substitution counting,
 * mate-overlap dedup, quality/conversion classification) lives here.  Python
 * only fills a cm_config and calls cm_run(), then reports the result.
 *
 * Built on the self-contained htslib subset from lh3/minipileup.
 */
#ifndef COUNTMUT_CORE_H
#define COUNTMUT_CORE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Modes */
#define CM_MODE_MUTATION 0
#define CM_MODE_BASE     1
#define CM_MODE_ALLELE   2

/* Engines */
#define CM_ENGINE_AUTO    0
#define CM_ENGINE_READWALK 1
#define CM_ENGINE_PILEUP   2

/* Strand processing */
#define CM_STRAND_BOTH    0
#define CM_STRAND_FORWARD 1
#define CM_STRAND_REVERSE 2

typedef struct {
    int32_t mode;            /* CM_MODE_* */
    int32_t engine;          /* CM_ENGINE_* */
    int32_t vcf;             /* allele mode: emit VCF */
    int     ref_base;        /* 'A'..'T' (mutation) */
    int     mut_base;
    int     ref_base2;       /* 0 = unset (alternative tagging) */
    int     mut_base2;
    int     pad;             /* motif half-window */
    int     save_rest;       /* emit o0/o1/o2 */
    int     min_mapq;
    int     min_baseq;
    int     max_sub;         /* NS cap, -1 = ignore */
    int     max_unc;         /* Zf cap, -1 = ignore */
    int     min_con;         /* Yf floor, -1 = ignore */
    int     trim_start;
    int     trim_end;
    int     min_allele_support;
    double  min_allele_frac;
    int     min_strand_support;
    int     min_depth;       /* base/allele: min site depth */
    int     mean_depth;
    int     count_indels;
    int     split_strand;
    int     strand_process;  /* CM_STRAND_* */
    int     max_depth;
    int     threads;
    int     flanking;        /* reference window w/o motif (pbr -k) */
    int     verbose;
    /* samtools-style read filtering */
    int     req_flags;       /* --rf/--incl-flags / reqflags */
    int     excl_flags;      /* --ff/--excl-flags / exclflags (default UNMAP|SEC|DUP|QCFAIL) */
    /* BED region restriction (pbr -b / -x) */
    const char *bedfile;     /* include-regions BED (or position list) */
    const char *exclude;     /* exclude-regions BED */
} cm_config;

/* Count `region` (NULL = whole file) of `bam` against `fa`, writing the output
 * in the configured mode.  Returns 0 on success, nonzero on error. */
int cm_run(const cm_config *cfg, const char *bam, const char *fa,
           const char *out_path, const char *region);

#ifdef __cplusplus
}
#endif

#endif /* COUNTMUT_CORE_H */
