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
        "  --mode M              mutation | base | allele (default mutation)\n"
        "  --engine E            read-walk | pileup | auto (default auto)\n"
        "  --vcf                 allele mode: emit VCF\n"
        "  --ref-base C          reference base (mutation)\n"
        "  --mut-base C          mutation base (mutation)\n"
        "  --pad N               motif half-window\n"
        "  --save-rest           emit o0/o1/o2\n"
        "  --min-mapq N --min-baseq N --max-sub N --max-unc N --min-con N\n"
        "  --trim-start N --trim-end N\n"
        "  --min-allele-support N --min-allele-frac F --min-strand-support N\n"
        "  --min-depth N --mean-depth N\n"
        "  --count-indels --split-strand\n"
        "  --strand S            both | forward | reverse\n"
        "  --max-depth N --threads N --flanking N\n");
}

static int mode_from(const char *s) {
    if (!strcmp(s, "base")) return CM_MODE_BASE;
    if (!strcmp(s, "allele")) return CM_MODE_ALLELE;
    return CM_MODE_MUTATION;
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
    cfg.mode = CM_MODE_MUTATION;
    cfg.engine = CM_ENGINE_AUTO;
    cfg.ref_base = 'A';
    cfg.mut_base = 'G';
    cfg.pad = 15;
    cfg.save_rest = 0;
    cfg.min_mapq = 0;
    cfg.min_baseq = 20;
    cfg.max_sub = 1;
    cfg.max_unc = 3;
    cfg.min_con = 1;
    cfg.trim_start = 2;
    cfg.trim_end = 2;
    cfg.min_allele_support = 1;
    cfg.min_allele_frac = 0.0;
    cfg.min_depth = 0;
    cfg.mean_depth = 0;
    cfg.count_indels = 0;
    cfg.split_strand = 1;
    cfg.strand_process = CM_STRAND_BOTH;
    cfg.max_depth = 8000;
    cfg.threads = 1;
    cfg.flanking = 0;
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
        {"ref-base", required_argument, 0, 1000},
        {"mut-base", required_argument, 0, 1001},
        {"pad", required_argument, 0, 'p'},
        {"save-rest", no_argument, 0, 1002},
        {"min-mapq", required_argument, 0, 'q'},
        {"min-baseq", required_argument, 0, 'Q'},
        {"max-sub", required_argument, 0, 1003},
        {"max-unc", required_argument, 0, 1004},
        {"min-con", required_argument, 0, 1005},
        {"trim-start", required_argument, 0, 1006},
        {"trim-end", required_argument, 0, 1007},
        {"min-allele-support", required_argument, 0, 1008},
        {"min-allele-frac", required_argument, 0, 1009},
        {"min-depth", required_argument, 0, 1010},
        {"mean-depth", required_argument, 0, 1011},
        {"count-indels", no_argument, 0, 1012},
        {"split-strand", no_argument, 0, 1013},
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
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0},
    };
    int c;
    while ((c = getopt_long(argc, argv, "b:f:o:r:m:e:p:q:Q:s:t:h", long_opts, NULL)) != -1) {
        switch (c) {
        case 2000: bam = optarg; break;
        case 'f': fa = optarg; break;
        case 'o': out = optarg; break;
        case 'r': region = optarg; break;
        case 'm': cfg.mode = mode_from(optarg); break;
        case 'e': cfg.engine = engine_from(optarg); break;
        case 'v': cfg.vcf = 1; break;
        case 'p': cfg.pad = atoi(optarg); break;
        case 'q': cfg.min_mapq = atoi(optarg); break;
        case 'Q': cfg.min_baseq = atoi(optarg); break;
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
        case 1000: cfg.ref_base = toupper((unsigned char)optarg[0]); break;
        case 1001: cfg.mut_base = toupper((unsigned char)optarg[0]); break;
        case 1002: cfg.save_rest = 1; break;
        case 1003: cfg.max_sub = atoi(optarg); break;
        case 1004: cfg.max_unc = atoi(optarg); break;
        case 1005: cfg.min_con = atoi(optarg); break;
        case 1006: cfg.trim_start = atoi(optarg); break;
        case 1007: cfg.trim_end = atoi(optarg); break;
        case 1008: cfg.min_allele_support = atoi(optarg); break;
        case 1009: cfg.min_allele_frac = atof(optarg); break;
        case 1010: cfg.min_depth = atoi(optarg); break;
        case 1011: cfg.mean_depth = atoi(optarg); break;
        case 1012: cfg.count_indels = 1; break;
        case 1013: cfg.split_strand = 1; break;
        case 1014: cfg.max_depth = atoi(optarg); break;
        case 1015: cfg.flanking = atoi(optarg); break;
        default: usage(); return 1;
        }
    }
    if (!bam || !fa) { usage(); return 1; }
    /* Engine resolution: auto -> pileup (single pass, general) */
    if (cfg.engine == CM_ENGINE_AUTO) cfg.engine = CM_ENGINE_PILEUP;
    if (cfg.engine == CM_ENGINE_PILEUP) {
        return cm_run(&cfg, bam, fa, out, region);
    }
    /* read-walk engine: fall back to pileup for now (correct output) */
    /* TODO: dedicated read-walk engine */
    fprintf(stderr, "[countmut_core] read-walk engine not yet available; using pileup\n");
    return cm_run(&cfg, bam, fa, out, region);
}
