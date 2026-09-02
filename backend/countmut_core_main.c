/* countmut_core_main.c -- CLI wrapper around cm_run() (the C countmut core).
 *
 * Compiled into the `countmut_core` binary.  Python shells out to this, so all
 * the computation stays in C and no htslib symbols leak into the Python process.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <getopt.h>
#include "countmut_core.h"
#include "countmut_expr.h"

/* samtools-style symbolic SAM flag name -> bit */
static int flag_bit(const char *name) {
    if (!strcasecmp(name, "PAIRED")) return 1;
    if (!strcasecmp(name, "PROPER_PAIR")) return 2;
    if (!strcasecmp(name, "UNMAP")) return 4;
    if (!strcasecmp(name, "MUNMAP")) return 8;
    if (!strcasecmp(name, "REVERSE")) return 16;
    if (!strcasecmp(name, "MREVERSE")) return 32;
    if (!strcasecmp(name, "READ1")) return 64;
    if (!strcasecmp(name, "READ2")) return 128;
    if (!strcasecmp(name, "SECONDARY")) return 256;
    if (!strcasecmp(name, "QCFAIL")) return 512;
    if (!strcasecmp(name, "DUP")) return 1024;
    if (!strcasecmp(name, "SUPPLEMENTARY")) return 2048;
    return 0;
}

/* parse a comma-separated flag spec that may combine symbolic names and hex/dec ints */
static int parse_flags(const char *spec) {
    int val = 0;
    char buf[512]; strncpy(buf, spec, 511); buf[511] = 0;
    char *tok = strtok(buf, ",");
    while (tok) {
        while (*tok == ' ') ++tok;
        if (tok[0] == '0' && (tok[1] == 'x' || tok[1] == 'X')) val |= (int)strtol(tok, NULL, 16);
        else {
            int b = flag_bit(tok);
            if (b) val |= b;
            else if (isdigit((unsigned char)tok[0])) val |= (int)strtol(tok, NULL, 10);
        }
        tok = strtok(NULL, ",");
    }
    return val;
}

/* parse samtools --input-fmt-option (comma-separated key[=value]) */
static void parse_input_fmt(const char *spec, int *req, int *excl, int *min_mq) {
    char buf[1024]; strncpy(buf, spec, 1023); buf[1023] = 0;
    char *tok = strtok(buf, ",");
    while (tok) {
        char *eq = strchr(tok, '=');
        if (!eq) { *req |= flag_bit(tok); tok = strtok(NULL, ","); continue; }
        *eq = 0; char *key = tok, *val = eq + 1;
        if (!strcmp(key, "reqflags")) *req = parse_flags(val);
        else if (!strcmp(key, "exclflags")) *excl = parse_flags(val);
        else if (!strcmp(key, "min-MQ") || !strcmp(key, "min-mapq")) *min_mq = atoi(val);
        tok = strtok(NULL, ",");
    }
}

static void usage(void) {
    fprintf(stderr,
        "Usage: countmut_core --bam FILE --fa FILE --out FILE [options]\n"
        "  --bam FILE            input BAM (required)\n"
        "  --fa FILE             reference FASTA (required)\n"
        "  --out FILE            output (or '-'/omit for stdout)\n"
        "  --region S            chr:start-end region (default: whole file)\n"
        "  --output-expr STR      -o output-row template: text + {expr} cells\n"
        "  --fmt-header STR       header line for a custom -o template\n"
        "  --engine E            read-walk | pileup | auto (default auto)\n"
        "  --vcf                 allele output: emit VCF\n"
        "  --min-mapq N --max-sub N --max-unc N --min-con N\n"
        "  --trim-fragment-start N --trim-fragment-end N   fragment 5'/3' trim\n"
        "  --trim-r1-end N --trim-r2-start N               read R1 3'-end / R2 5'-start trim\n"
        "  --min-allele-support N --min-allele-frac F --min-strand-support N\n"
        "  --min-depth N --mean-depth N\n"
        "  --motif-pad N         {motif} reference window: 2*N+1 bases (0 = the base only)\n"
        "  --count-indels [--strandless]\n"
        "  --strand S            both | forward | reverse\n"
        "  --read-expr EXPR       -e Lua read filter (evaluated per base)\n"
        "  --pile-expr EXPR       -p Lua site filter (evaluated per site)\n"
        "  --verbose               real-time per-region progress on stderr\n"
        "  --max-depth N          pileup per-position depth cap (0 = unlimited)\n"
        "  --threads N --flanking N\n");
}

static int engine_from(const char *s) {
    if (!strcmp(s, "read-walk")) return CM_ENGINE_READWALK;
    if (!strcmp(s, "pileup")) return CM_ENGINE_PILEUP;
    return CM_ENGINE_AUTO;
}
static int strand_from(const char *s) {
    if (!strcmp(s, "forward")) return CM_STRAND_FORWARD;
    if (!strcmp(s, "reverse")) return CM_STRAND_REVERSE;
    return CM_STRAND_BOTH;
}

int main(int argc, char **argv) {
    cm_config cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.out = CM_OUT_COMPOSITION;
    cfg.engine = CM_ENGINE_AUTO;
    cfg.min_mapq = 0;
    cfg.min_baseq = 0;     /* quality QC is a -e filter (bq >= N): all counted bases are tier x2 */
    cfg.max_sub = -1;      /* read filters + trimming live in -e, not defaults */
    cfg.max_unc = -1;
    cfg.min_con = -1;
    cfg.trim_fragment_start = 0;
    cfg.trim_fragment_end = 0;
    cfg.trim_r1_end = 0;
    cfg.trim_r2_start = 0;
    cfg.min_allele_support = 1;
    cfg.min_allele_frac = 0.0;
    cfg.min_depth = 0;
    cfg.mean_depth = 0;
    cfg.count_indels = 0;
    cfg.strandless = 0;      /* base/allele: per-strand '+'/'-' by default; --strandless collapses */
    cfg.strand_process = CM_STRAND_BOTH;
    cfg.max_depth = 0;       /* 0 = unlimited (count all reads) */
    cfg.threads = 1;
    cfg.flanking = 0;
    cfg.pad = 0;              /* {motif} window: 2*pad+1 ref bases (0 = the base only) */
    cfg.req_flags = 0;
    cfg.excl_flags = 1796; /* samtools default: UNMAP|SECONDARY|QCFAIL|DUP */
    cfg.bedfile = NULL;
    cfg.exclude = NULL;

    const char *bam = NULL, *fa = NULL, *out = NULL, *region = NULL;
    static struct option long_opts[] = {
        {"bam", required_argument, 0, 2000},
        {"fa", required_argument, 0, 'f'},
        {"out", required_argument, 0, 'o'},
        {"region", required_argument, 0, 'r'},
        {"mode", required_argument, 0, 'm'},
        {"engine", required_argument, 0, 'e'},
        {"vcf", no_argument, 0, 'v'},
        {"min-mapq", required_argument, 0, 'q'},
        {"max-sub", required_argument, 0, 1003},
        {"max-unc", required_argument, 0, 1004},
        {"min-con", required_argument, 0, 1005},
        {"trim-fragment-start", required_argument, 0, 1006},
        {"trim-fragment-end", required_argument, 0, 1007},
        {"trim-r1-end", required_argument, 0, 1017},
        {"trim-r2-start", required_argument, 0, 1018},
        {"min-allele-support", required_argument, 0, 1008},
        {"min-allele-frac", required_argument, 0, 1009},
        {"min-depth", required_argument, 0, 1010},
        {"mean-depth", required_argument, 0, 1011},
        {"count-indels", no_argument, 0, 1012},
        {"strandless", no_argument, 0, 1013},
        {"strand", required_argument, 0, 's'},
        {"max-depth", required_argument, 0, 1014},
        {"threads", required_argument, 0, 't'},
        {"flanking", required_argument, 0, 'k'},
        {"bedfile", required_argument, 0, 'b'},
        {"exclude", required_argument, 0, 'x'},
        {"incl-flags", required_argument, 0, 1100},
        {"rf", required_argument, 0, 1100},
        {"excl-flags", required_argument, 0, 1101},
        {"ff", required_argument, 0, 1101},
        {"input-fmt-option", required_argument, 0, 1102},
        {"mate-fix", no_argument, 0, 1016},
        {"motif-pad", required_argument, 0, 1019},
        {"read-expr", required_argument, 0, 2002},
        {"pile-expr", required_argument, 0, 2003},
        {"output-expr", required_argument, 0, 2004},
        {"fmt-header", required_argument, 0, 2005},
        {"verbose", no_argument, 0, 'V'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0},
    };
    int c;
    while ((c = getopt_long(argc, argv, "b:f:o:r:e:q:s:t:h", long_opts, NULL)) != -1) {
        switch (c) {
        case 2000: bam = optarg; break;
        case 2002: cfg.read_expr = optarg; break;
        case 2003: cfg.pile_expr = optarg; break;
        case 2004: cfg.output_expr = optarg; break;
        case 2005: cfg.fmt_header = optarg; break;
        case 'V': cfg.verbose = 1; break;
        case 'f': fa = optarg; break;
        case 'o': out = optarg; break;
        case 'r': region = optarg; break;
        case 'e': cfg.engine = engine_from(optarg); break;
        case 'v': cfg.vcf = 1; break;
        case 'q': cfg.min_mapq = atoi(optarg); break;
        case 's': cfg.strand_process = strand_from(optarg); break;
        case 't': cfg.threads = atoi(optarg); break;
        case 'b': cfg.bedfile = optarg; break;
        case 'x': cfg.exclude = optarg; break;
        case 'k': cfg.flanking = atoi(optarg); break;
        case 'h': usage(); return 0;
        case 1100: cfg.req_flags |= parse_flags(optarg); break;
        case 1101: cfg.excl_flags = parse_flags(optarg); break;
        case 1102: parse_input_fmt(optarg, &cfg.req_flags, &cfg.excl_flags, &cfg.min_mapq); break;
        case 1016: /* --mate-fix: overlap dedup is always enabled for correctness */ break;
        case 1003: cfg.max_sub = atoi(optarg); break;
        case 1004: cfg.max_unc = atoi(optarg); break;
        case 1005: cfg.min_con = atoi(optarg); break;
        case 1006: cfg.trim_fragment_start = atoi(optarg); break;
        case 1007: cfg.trim_fragment_end = atoi(optarg); break;
        case 1017: cfg.trim_r1_end = atoi(optarg); break;
        case 1018: cfg.trim_r2_start = atoi(optarg); break;
        case 1008: cfg.min_allele_support = atoi(optarg); break;
        case 1009: cfg.min_allele_frac = atof(optarg); break;
        case 1010: cfg.min_depth = atoi(optarg); break;
        case 1011: cfg.mean_depth = atoi(optarg); break;
        case 1012: cfg.count_indels = 1; break;
        case 1013: cfg.strandless = 1; break;
        case 1014: cfg.max_depth = atoi(optarg); break;
        case 1015: cfg.flanking = atoi(optarg); break;
        case 1019: cfg.pad = atoi(optarg); break;
        default: usage(); return 1;
        }
    }
    if (!bam || !fa) { usage(); return 1; }
    /* Validate -e / -p / -o expressions up-front (syntax error -> exit 2). */
    if (!cm_expr_valid(cfg.read_expr, cfg.pile_expr, cfg.output_expr)) return 2;
    /* Output format: --vcf -> allele (VCF) output; otherwise composition.
     * There is no `--mode`; everything else is a --output-format template. */
    cfg.out = cfg.vcf ? CM_OUT_ALLELE : CM_OUT_COMPOSITION;
    /* Engine: pileup for the per-position counting (assume no `--engine` given
     * or auto); both engines emit identical rows. */
    if (cfg.engine == CM_ENGINE_AUTO) cfg.engine = CM_ENGINE_PILEUP;
    return cm_run(&cfg, bam, fa, out, region);
}
