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

/* Engines */
#define CM_ENGINE_AUTO    0
#define CM_ENGINE_READWALK 1
#define CM_ENGINE_PILEUP   2

/* -e read-expression status categories.  Slot 0 is the DEFAULT category
 * (used when no -e filter is set, or the filter returns true); slots 1..K
 * are the user-declared --status names (in declaration order).  site_t keeps
 * one [strand][category][base] count matrix, so the -p / output expressions
 * can split the per-site counts by category (e.g. legacy u0/u1/u2 columns).
 * Bumped from 3 -> 4 to fit the default + the three legacy tiers
 * (low quality / high conversion / insufficient conversion). */
#define CM_CAT_MAX 4
#define CM_CAT_NAME_MAX 16

/* per-site accumulator: counts per biological strand, per -e status
 * category, per base (0=A,1=C,2=G,3=T,4=N) */
typedef struct {
    int cnt[2][CM_CAT_MAX][5]; /* [strand][category][base] */
    int ins[2], del[2], refskip[2], fail[2];
} site_t;

/* Strand processing */
#define CM_STRAND_BOTH    0
#define CM_STRAND_FORWARD 1
#define CM_STRAND_REVERSE 2

#define CM_OUT_CONVERSION  0  /* legacy/unused: u/m conversion view (dead)   */
#define CM_OUT_COMPOSITION 1  /* output format: per-base composition (default) */
#define CM_OUT_ALLELE      2  /* output format: ref/alt / VCF                */

typedef struct {
    int32_t out;             /* CM_OUT_* (output format) */
    int32_t engine;          /* CM_ENGINE_* */
    int32_t vcf;             /* allele output: emit VCF */
    int     ref_base;        /* legacy/unused: always 0 (no --ref-base flag) */
    int     mut_base;
    int     ref_base2;       /* legacy/unused: always 0 */
    int     mut_base2;
    int     pad;             /* {motif} reference window: 2*pad+1 bases (--motif-pad) */
    int     save_rest;       /* legacy/unused */
    const char *output_expr; /* -o output-row template (overrides the built-in format) */
    const char *fmt_header;  /* header line for a custom output template ("" = none) */
    int     min_mapq;
    int     min_baseq;
    int     max_sub;         /* NS cap, -1 = ignore */
    int     max_unc;         /* Zf cap, -1 = ignore */
    int     min_con;         /* Yf floor, -1 = ignore */
    int     trim_fragment_start;  /* fragment 5' end trim */
    int     trim_fragment_end;    /* fragment 3' end trim */
    int     trim_r1_end;          /* read1 (R1) 3'-query-end trim */
    int     trim_r2_start;        /* read2 (R2) 5'-query-start trim */
    int     min_allele_support;
    double  min_allele_frac;
    int     min_strand_support;
    int     min_depth;       /* base/allele: min site depth */
    int     mean_depth;
    int     count_indels;
    int     strandless;      /* 1 = collapse +/- strands (base/allele); 0 = per-strand (default) */
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
    /* Lua filter expressions (-e read-level, -p site-level) */
    const char *read_expr;
    const char *pile_expr;
} cm_config;

/* Count `region` (NULL = whole file) of `bam` against `fa`, writing the output
 * in the configured mode.  Returns 0 on success, nonzero on error. */
int cm_run(const cm_config *cfg, const char *bam, const char *fa,
           const char *out_path, const char *region);

#ifdef __cplusplus
}
#endif

#endif /* COUNTMUT_CORE_H */
