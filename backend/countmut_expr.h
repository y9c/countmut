/* countmut_expr.h -- Lua filter expressions for the C countmut core.
 *
 * Evaluates the -e (read-level, once per aligned base) and -p (site-level,
 * once per reported site) filters inside the C core, following pbr's approach
 * of embedding Lua.  When neither filter is set the evaluator is NULL and the
 * hot paths skip it entirely (zero overhead).
 */
#ifndef COUNTMUT_EXPR_H
#define COUNTMUT_EXPR_H

#include "sam.h"

typedef struct cm_expr cm_expr;

/* Create an evaluator for the given -e / -p Lua expressions (either may be
 * NULL).  Returns NULL when there is nothing to evaluate.  On a compile error
 * the offending message is printed to stderr and that filter is disabled;
 * use cm_expr_valid() first to turn syntax errors into a hard exit. */
cm_expr *cm_expr_new(const char *read_expr, const char *pile_expr);

/* Validate both expressions (compile-only).  Returns 1 if both valid or both
 * absent (or NULL -> 1); 0 if a compile error was reported to stderr. */
int cm_expr_valid(const char *read_expr, const char *pile_expr);

void cm_expr_free(cm_expr *x);

/* True if a read-level (-e) filter is configured and compiled. */
int cm_expr_has_read(const cm_expr *x);
/* True if the -e filter depends only on read-level values (not per-base
 * qpos/bq/base/ref/dist) -- callers may then evaluate it once per read. */
int cm_expr_read_constant(const cm_expr *x);
/* True if a site-level (-p) filter is configured and compiled. */
int cm_expr_has_pile(const cm_expr *x);

/* Evaluate the -e filter for one aligned base of read `b`.
 * `rname`/`mrname` are the contig names of the read and its mate (may be ""),
 * `qpos` the query position, `strand_sign` is +1/-1 (biological, paired-aware)
 * and `ref_base` the reference base the base aligns over.  Returns 1 = keep
 * the base, 0 = reject.  Always 1 if no read filter. */
int cm_expr_read(cm_expr *x, const bam1_t *b, const char *rname, const char *mrname,
                 int qpos, int strand_sign, char ref_base);

/* Evaluate the -p filter for one site.  cnt[5] = A/C/G/T/N totals across all
 * strands and quality tiers; ins/del/rs/fl are the indel/ref-skip/fail counts;
 * motif is the reference window string (may be NULL/empty outside mutation
 * mode).  Returns 1 = keep the site, 0 = omit it. */
int cm_expr_pile(cm_expr *x, int64_t pos, char ref_ch, const char *motif,
                 const int cnt[5], int ins, int del, int rs, int fl);

#endif /* COUNTMUT_EXPR_H */
