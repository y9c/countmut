/* countmut_expr.c -- Lua filter expressions (-e / -p) for the C countmut core.
 *
 * Follows pbr: embed Lua 5.4, compile each filter string as a Lua chunk once,
 * and evaluate it with the read / site values exposed as global variables (the
 * countmut flat namespace: mapq, bq, flags, strand, qname, pos, dist5/dist3,
 * tag('XX')...) and as the `read` / `pile` tables (the pbr object style).
 * A filter that evaluates to false rejects the base (read) or omits the site
 * (pile).  Runtime errors are treated as false (never propagate).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <regex.h>

#include "countmut_expr.h"
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#define CUR_READ_KEY "countmut.request_read"

struct cm_expr {
    lua_State *L;
    int read_ref;   /* LUA_NOREF if no -e */
    int pile_ref;   /* LUA_NOREF if no -p */
    int need_seq;    /* read seq / read.sequence */
    int need_qual;   /* raw quality string */
    int need_library;/* library / read.library / LB */
    int need_qname, need_rname, need_mrname;
    int need_mapq, need_flag, need_strand, need_pos, need_endpos, need_tid;
    int need_mate, need_tlen, need_qlen, need_rlen, need_cigar, need_flags2;
    int need_qpos, need_bq, need_base, need_ref, need_dist;  /* per-base */
    int need_aux;    /* expression calls tag()/exists()/n5()/n3() (needs CUR_READ_KEY) */
    int need_read_tbl, need_pile_tbl;  /* expression uses `read.x` / `pile.x` dotted access */
    int read_reg, pile_reg;            /* luaL_ref of the read / pile tables */
};

#if LUA_VERSION_NUM < 502
#error "countmut_expr needs Lua >= 5.2"
#endif

/* ---- aux tag values (bam_aux_get returns [type][payload...]) ---- */
static void push_aux_value(lua_State *L, const uint8_t *s) {
    switch (*s) {
    case 'A': case 'a':
        lua_pushinteger(L, bam_aux2A(s)); break;
    case 'c': case 'C': case 's': case 'S': case 'i': case 'I':
        lua_pushinteger(L, bam_aux2i(s)); break;
    case 'f': case 'F': case 'd': case 'D':
        lua_pushnumber(L, bam_aux2f(s)); break;
    case 'Z': case 'z': case 'H': case 'h':
        lua_pushstring(L, (const char *)(s + 1)); break;
    default:
        lua_pushnil(L); break;
    }
}

static const bam1_t *current_read(lua_State *L) {
    const bam1_t *b;
    lua_getfield(L, LUA_REGISTRYINDEX, CUR_READ_KEY);
    b = (const bam1_t *)lua_touserdata(L, -1);
    lua_pop(L, 1);
    return b;
}

/* tag('XX'): read an aux tag from the current read (nil if missing/unreadable) */
static int l_tag(lua_State *L) {
    const char *name = luaL_checkstring(L, 1);
    const bam1_t *b = current_read(L);
    if (b == NULL || name == NULL || name[0] == '\0' || name[2] != '\0') {
        lua_pushnil(L);
        return 1;
    }
    uint8_t *aux = bam_aux_get(b, name);
    if (aux == NULL) { lua_pushnil(L); return 1; }
    push_aux_value(L, aux);
    return 1;
}

/* exists('XX'): is the aux tag present on the current read? */
static int l_exists(lua_State *L) {
    const char *name = luaL_checkstring(L, 1);
    const bam1_t *b = current_read(L);
    int ok = (b != NULL && name != NULL && name[0] != '\0' && name[2] == '\0'
              && bam_aux_get(b, name) != NULL);
    lua_pushboolean(L, ok);
    return 1;
}

/* ---- string / regex helpers (samtools FILTER EXPRESSIONS style) ---- */
static int l_re_match(lua_State *L) {
    const char *subj = luaL_checkstring(L, 1);
    const char *pat = luaL_checkstring(L, 2);
    regex_t re; int matched = 0;
    if (subj && pat && regcomp(&re, pat, REG_EXTENDED | REG_NOSUB) == 0) {
        matched = (regexec(&re, subj, 0, NULL, 0) == 0);
        regfree(&re);
    }
    lua_pushboolean(L, matched);
    return 1;
}
static int l_slen(lua_State *L) {
    size_t n;
    luaL_checklstring(L, 1, &n);
    lua_pushinteger(L, (lua_Integer)n);
    return 1;
}
static int l_smin(lua_State *L) {
    size_t n; const char *s = luaL_checklstring(L, 1, &n);
    int m = -1;
    for (size_t i = 0; i < n; ++i) { int v = (unsigned char)s[i]; if (i == 0 || v < m) m = v; }
    lua_pushinteger(L, m);
    return 1;
}
static int l_smax(lua_State *L) {
    size_t n; const char *s = luaL_checklstring(L, 1, &n);
    int m = -1;
    for (size_t i = 0; i < n; ++i) { int v = (unsigned char)s[i]; if (i == 0 || v > m) m = v; }
    lua_pushinteger(L, m);
    return 1;
}
static int l_savg(lua_State *L) {
    size_t n; const char *s = luaL_checklstring(L, 1, &n);
    if (n == 0) { lua_pushnumber(L, 0.0 / 0.0); return 1; }   /* NaN like samtools */
    double sum = 0;
    for (size_t i = 0; i < n; ++i) sum += (unsigned char)s[i];
    lua_pushnumber(L, sum / (double)n);
    return 1;
}
static const char *read_seq_str(const bam1_t *b, char *buf, size_t cap) {
    if (!b || b->core.l_qseq == 0) { buf[0] = 0; return buf; }
    size_t l = (size_t)b->core.l_qseq;
    if (l + 1 > cap) l = cap - 1;
    static const char ntc[] = "=ACMGRSVTWYHKDBN";
    for (size_t i = 0; i < l; ++i) {
        int v = (int)bam_seqi(bam_get_seq(b), (int)i);
        buf[i] = (v < 16 && ntc[v]) ? ntc[v] : 'N';
    }
    buf[l] = 0;
    return buf;
}
/* N-proportion in the first(n5)/last(n3) n bases of the current read's sequence */
static int l_n_prop(lua_State *L, int from_end) {
    int n = (int)luaL_checkinteger(L, 1);
    const bam1_t *b = current_read(L);
    if (!b) { lua_pushnumber(L, 0); return 1; }
    int lq = (int)b->core.l_qseq;
    if (lq == 0) { lua_pushnumber(L, 0); return 1; }
    if (n < 0) n = 0;
    if (n > lq) n = lq;
    int ns = 0;
    for (int i = 0; i < n; ++i) {
        int qi = from_end ? (lq - n + i) : i;
        int v = (int)bam_seqi(bam_get_seq(b), qi);
        if (v == 15) ns++;   /* N */
    }
    lua_pushnumber(L, (double)ns / (double)lq);
    return 1;
}
static int l_n5(lua_State *L) { return l_n_prop(L, 0); }
static int l_n3(lua_State *L) { return l_n_prop(L, 1); }

static int is_rev(const bam1_t *b) { return (b->core.flag & 16) != 0; }

/* htslib-style math helpers (sqrt/log/exp/pow) -- thin wrappers over math.* */
static int l_f1(lua_State *L, const char *f) {
    luaL_checknumber(L, 1);
    lua_getglobal(L, "math");
    lua_getfield(L, -1, f);
    lua_pushvalue(L, 1);
    lua_call(L, 1, 1);
    return 1;
}
static int l_sqrt(lua_State *L) { return l_f1(L, "sqrt"); }
static int l_log (lua_State *L) { return l_f1(L, "log"); }
static int l_exp (lua_State *L) { return l_f1(L, "exp"); }
static int l_pow (lua_State *L) {
    luaL_checknumber(L, 1);
    luaL_checknumber(L, 2);
    lua_getglobal(L, "math");
    lua_getfield(L, -1, "pow");
    lua_pushvalue(L, 1);
    lua_pushvalue(L, 2);
    lua_call(L, 2, 1);
    return 1;
}

static void ref_span(const bam1_t *b, int *rlen) {
    int r = 0;
    const uint32_t *cig = bam_get_cigar(b);
    for (int i = 0; i < b->core.n_cigar; ++i) {
        int op = (int)bam_cigar_op(cig[i]);
        if (op == 0 || op == 2 || op == 3 || op == 7 || op == 8)
            r += (int)bam_cigar_oplen(cig[i]);
    }
    *rlen = r;
}

/* cigar-derived read stats (samtools sclen/hclen/qlen/rlen + indel & soft-clip
 * ends, pbr soft_clips_5_prime/3_prime / indel_count) */
static void cigar_stats(const bam1_t *b, int *sclen, int *hclen, int *n_indel,
                        int *soft5, int *soft3) {
    int sc = 0, hc = 0, ni = 0, s5 = 0, s3 = 0;
    const uint32_t *cig = bam_get_cigar(b);
    for (int i = 0; i < b->core.n_cigar; ++i) {
        int op = (int)bam_cigar_op(cig[i]);
        int len = (int)bam_cigar_oplen(cig[i]);
        switch (op) {
        case 4: sc += len; break;
        case 5: hc += len; break;
        case 1: case 2: ni += 1; break;
        default: break;
        }
    }
    if (b->core.n_cigar >= 1 && (int)bam_cigar_op(cig[0]) == 4) s5 = (int)bam_cigar_oplen(cig[0]);
    if (b->core.n_cigar >= 2 && (int)bam_cigar_op(cig[b->core.n_cigar - 1]) == 4)
        s3 = (int)bam_cigar_oplen(cig[b->core.n_cigar - 1]);
    *sclen = sc; *hclen = hc; *n_indel = ni; *soft5 = s5; *soft3 = s3;
}

/* function globals that must never be overwritten by a field with the same name
 * (e.g. `length` is both read.length and the htslib-style length() helper) */
static int is_reserved(const char *name) {
    static const char *fns[] = {"tag","exists","re_match","re_find","slen","smin",
                                "smax","savg","length","min","max","avg","n5","n3"};
    for (size_t i = 0; i < sizeof(fns) / sizeof(fns[0]); ++i)
        if (strcmp(name, fns[i]) == 0) return 1;
    return 0;
}

/* set one scalar into the global namespace, and INTO the `read`/`pile` table
 * only when the expression actually uses dotted access (`read.x` / `pile.x`).
 * tbl_reg = luaL_ref of the table (cheaper than a string-keyed getglobal per
 * field), set_tbl = x->need_{read,pile}_tbl.  Most bare predicates (`bq >= 20`)
 * therefore skip the table-population entirely. */
static void put_int(lua_State *L, int tbl_reg, int set_tbl, const char *name, lua_Integer v) {
    if (!is_reserved(name)) { lua_pushinteger(L, v); lua_setglobal(L, name); }
    if (set_tbl) {
        lua_rawgeti(L, LUA_REGISTRYINDEX, tbl_reg);
        lua_pushinteger(L, v); lua_setfield(L, -2, name);
        lua_pop(L, 1);
    }
}
static void put_str(lua_State *L, int tbl_reg, int set_tbl, const char *name, const char *v) {
    if (!is_reserved(name)) { lua_pushstring(L, v); lua_setglobal(L, name); }
    if (set_tbl) {
        lua_rawgeti(L, LUA_REGISTRYINDEX, tbl_reg);
        lua_pushstring(L, v); lua_setfield(L, -2, name);
        lua_pop(L, 1);
    }
}
/* read/pile-scoped wrappers used from the per-base / per-site setters */
static void r_int(cm_expr *x, const char *name, lua_Integer v) { put_int(x->L, x->read_reg, x->need_read_tbl, name, v); }
static void r_str(cm_expr *x, const char *name, const char *v) { put_str(x->L, x->read_reg, x->need_read_tbl, name, v); }
static void p_int(cm_expr *x, const char *name, lua_Integer v) { put_int(x->L, x->pile_reg, x->need_pile_tbl, name, v); }
static void p_str(cm_expr *x, const char *name, const char *v) { put_str(x->L, x->pile_reg, x->need_pile_tbl, name, v); }

static int run_chunk(lua_State *L, int ref) {
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    if (lua_pcall(L, 0, 1, 0) != LUA_OK) {
        fprintf(stderr, "[countmut] expression error: %s\n", lua_tostring(L, -1));
        lua_pop(L, 1);
        return 0;                       /* reject on error, never propagate */
    }
    int keep = lua_toboolean(L, -1);
    lua_pop(L, 1);
    return keep;
}

/* Rewrite `IDENT =~ "re"` / `IDENT !~ "re"` (samtools regex grammar, IDENT a
 * dotted name such as `rname`, `read.qname` or `tag('RG')` is NOT auto-rewritten
 * -- use re_match(field, "pat") for those) into re_match(IDENT, "re") /
 * not re_match(IDENT, "re").  Non-matching patterns are left untouched. */
static char *translate_regex(const char *src) {
    size_t n = strlen(src);
    char *out = (char *)malloc(2 * n + 64);
    if (out == NULL) return NULL;
    size_t j = 0, i = 0;
    while (i < n) {
        char c = src[i];
        if (isalpha((unsigned char)c) || c == '_') {
            size_t s = i;
            while (i < n && (isalnum((unsigned char)src[i]) || src[i] == '_' || src[i] == '.'))
                i++;
            /* identifier = src[s..i) */
            size_t k = i;
            while (k < n && (src[k] == ' ' || src[k] == '\t')) k++;
            int neg = 0;
            if (k + 1 < n && src[k] == '=' && src[k + 1] == '~') neg = 0;
            else if (k + 1 < n && src[k] == '!' && src[k + 1] == '~') neg = 1;
            else {
                memcpy(out + j, src + s, i - s); j += i - s;
                continue;
            }
            size_t p = k + 2;
            while (p < n && (src[p] == ' ' || src[p] == '\t')) p++;
            if (p >= n || (src[p] != '"' && src[p] != '\'')) {
                memcpy(out + j, src + s, i - s); j += i - s;
                continue;   /* complex rhs; leave as-is */
            }
            char q = src[p];
            size_t qs = p, qe = p;
            for (;;) {
                if (qe >= n) { qe = p + 1; break; }
                if (src[qe] == q && qe > p) { qe++; break; }
                if (src[qe] == '\\' && qe + 1 < n) qe += 2; else qe++;
            }
            if (neg) { memcpy(out + j, "not re_match(", 13); j += 13; }
            else { memcpy(out + j, "re_match(", 9); j += 9; }
            memcpy(out + j, src + s, i - s); j += i - s;
            out[j++] = ','; out[j++] = ' ';
            memcpy(out + j, src + qs, qe - qs); j += qe - qs;
            out[j++] = ')';
            i = qe;
            continue;
        }
        out[j++] = c;
        i++;
    }
    out[j] = 0;
    return out;
}

/* countmut/Python/samtools-style expressions use !=, &&, ||, !, [XX] tags and
 * flag.NAME; Lua uses ~=, and, or, not, tag('XX'), flags & BIT.  Translate
 * (outside string literals) so samtools/pbr-style -e / -p expressions work. */
static char *translate_ops(const char *src0) {
    char *src = translate_regex(src0);
    if (src == NULL) return NULL;
    size_t n = strlen(src);
    char *out = (char *)malloc(2 * n + 64);
    if (out == NULL) { free(src); return NULL; }
    size_t j = 0;
    char quote = 0;
    for (size_t i = 0; i < n; ++i) {
        char c = src[i];
        if (quote) {
            out[j++] = c;
            if (c == '\\' && i + 1 < n) { out[j++] = src[++i]; }
            else if (c == quote) quote = 0;
            continue;
        }
        if (c == '\'' || c == '"') { quote = c; out[j++] = c; continue; }
        if (c == '!' && i + 1 < n && src[i + 1] == '=') { out[j++] = '~'; out[j++] = '='; i += 1; continue; }
        if (c == '&' && i + 1 < n && src[i + 1] == '&') { out[j++] = 'a'; out[j++] = 'n'; out[j++] = 'd'; i += 1; continue; }
        if (c == '|' && i + 1 < n && src[i + 1] == '|') { out[j++] = 'o'; out[j++] = 'r'; i += 1; continue; }
        if (c == '!') { out[j++] = 'n'; out[j++] = 'o'; out[j++] = 't'; out[j++] = ' '; continue; }
        if (c == '[' && i + 3 < n && isalnum((unsigned char)src[i + 1])
            && isalnum((unsigned char)src[i + 2]) && src[i + 3] == ']') {
            /* inside exists( ... ) -> exists('XX') (the htslib idiom);
             * anywhere else -> tag('XX') */
            if (j >= 7 && memcmp(out + j - 7, "exists(", 7) == 0) {
                out[j++] = '\''; out[j++] = src[i + 1]; out[j++] = src[i + 2];
                out[j++] = '\'';
            } else {
                out[j++] = 't'; out[j++] = 'a'; out[j++] = 'g'; out[j++] = '(';
                out[j++] = '\''; out[j++] = src[i + 1]; out[j++] = src[i + 2];
                out[j++] = '\''; out[j++] = ')';
            }
            i += 3;
            continue;
        }
        if (c == 'f' && n - i > 5 && memcmp(src + i, "flag.", 5) == 0) {
            static const struct { const char *n; int bit; } fl[] = {
                {"paired",1},{"proper_pair",2},{"unmap",4},{"munmap",8},
                {"reverse",16},{"mreverse",32},{"read1",64},{"read2",128},
                {"secondary",256},{"qcfail",512},{"dup",1024},{"supplementary",2048},
            };
            const char *nm = src + i + 5;
            size_t nl = 0; while (nm[nl] && (isalnum((unsigned char)nm[nl]) || nm[nl] == '_')) nl++;
            int bit = 0;
            for (size_t f = 0; f < sizeof(fl) / sizeof(fl[0]); ++f)
                if (nl == strlen(fl[f].n) && memcmp(nm, fl[f].n, nl) == 0) { bit = fl[f].bit; break; }
            if (bit) {
                out[j++] = 'f'; out[j++] = 'l'; out[j++] = 'a'; out[j++] = 'g';
                out[j++] = 's'; out[j++] = ' '; out[j++] = '&'; out[j++] = ' ';
                size_t b = 0;
                char bb[16]; int tmp = bit;
                if (tmp == 0) bb[b++] = '0';
                while (tmp) { bb[b++] = (char)('0' + tmp % 10); tmp /= 10; }
                while (b > 0) out[j++] = bb[--b];
                i += 5 + nl - 1;
                continue;
            }
        }
        out[j++] = c;
    }
    out[j] = 0;
    free(src);
    return out;
}

/* prefix-tolerant name search: `flag` matches `flags`, `read.qlen` matches
 * `qlen`; matching too much is safe (an extra field is set), matching too
 * little silently leaves a field nil -- so we are deliberately liberal (a
 * token after `.` or `_` is still a field / name we should populate). */
static int has_token(const char *s, const char *tok) {
    size_t tl = strlen(tok);
    for (const char *p = s; (p = strstr(p, tok)); p += tl) {
        if (p == s || !(isalnum((unsigned char)p[-1]) || p[-1] == '_'))
            return 1;
    }
    return 0;
}

static int compile_chunk(lua_State *L, const char *expr, const char *what,
                         int *need_seq, int *need_qual, int *need_lib,
                         int *need_qname, int *need_rname, int *need_mrname,
                         int *need_mapq, int *need_flag, int *need_strand,
                         int *need_pos, int *need_endpos, int *need_tid,
                         int *need_mate, int *need_tlen, int *need_qlen,
                         int *need_rlen, int *need_cigar, int *need_flags2,
                         int *need_qpos, int *need_bq, int *need_base,
                         int *need_ref, int *need_dist, int *need_aux,
                         int *need_read_tbl, int *need_pile_tbl) {
    if (expr == NULL || expr[0] == '\0') return LUA_NOREF;
    char *t = translate_ops(expr);
    const char *src = t ? t : expr;
    if (need_seq)     *need_seq     = has_token(src, "sequence") || has_token(src, "seq");
    if (need_qual)    *need_qual    = has_token(src, "qual");
    if (need_lib)     *need_lib     = has_token(src, "library") || has_token(src, "LB");
    if (need_qname)   *need_qname   = has_token(src, "qname");
    if (need_rname)   *need_rname   = has_token(src, "rname");
    if (need_mrname)  *need_mrname  = has_token(src, "mrname") || has_token(src, "rnext");
    if (need_mapq)    *need_mapq    = has_token(src, "mapq") || has_token(src, "MAPQ")
                                      || has_token(src, "mapping_quality");
    if (need_flag)    *need_flag    = has_token(src, "flag");
    if (need_strand)  *need_strand  = has_token(src, "strand") || has_token(src, "STRAND");
    if (need_pos)     *need_pos     = has_token(src, "pos") || has_token(src, "POS") || has_token(src, "start");
    if (need_endpos)  *need_endpos  = has_token(src, "endpos") || has_token(src, "stop");
    if (need_tid)     *need_tid     = has_token(src, "tid") || has_token(src, "refid");
    if (need_mate)    *need_mate    = has_token(src, "mtid") || has_token(src, "mpos") || has_token(src, "pnext");
    if (need_tlen)    *need_tlen    = has_token(src, "tlen") || has_token(src, "insert_size");
    if (need_qlen)    *need_qlen    = has_token(src, "qlen") || has_token(src, "length") || has_token(src, "LEN");
    if (need_rlen)    *need_rlen    = has_token(src, "rlen");
    if (need_cigar)   *need_cigar   = has_token(src, "sclen") || has_token(src, "hclen")
                                      || has_token(src, "indel_count") || has_token(src, "n_indel")
                                      || has_token(src, "soft_clips") || has_token(src, "ncigar");
    if (need_flags2)  *need_flags2  = has_token(src, "is_reverse") || has_token(src, "is_paired")
                                      || has_token(src, "r1") || has_token(src, "r2");
    if (need_qpos)    *need_qpos    = has_token(src, "qpos") || has_token(src, "QPOS");
    if (need_bq)      *need_bq      = has_token(src, "bq") || has_token(src, "BQ") || has_token(src, "baseq");
    if (need_base)    *need_base    = has_token(src, "base");
    if (need_ref)     *need_ref     = has_token(src, "ref");
    if (need_dist)    *need_dist    = has_token(src, "dist5") || has_token(src, "dist3")
                                      || has_token(src, "DIST5") || has_token(src, "DIST3")
                                      || has_token(src, "distance_from_5prime")
                                      || has_token(src, "distance_from_3prime");
    if (need_aux)     *need_aux     = has_token(src, "tag(") || has_token(src, "exists(")
                                      || has_token(src, "n5(") || has_token(src, "n3(")
                                      || has_token(src, "n_proportion");
    if (need_read_tbl && strcmp(what, "read") == 0) *need_read_tbl = has_token(src, "read.");
    if (need_pile_tbl && strcmp(what, "pile") == 0) *need_pile_tbl = has_token(src, "pile.");
    /* Always load as `return (...)` so a bare predicate or a bare call (e.g.
     * `rname =~ 'x'` -> re_match(...)) yields a VALUE, not a nil-returning
     * call-statement.  A source that already begins with `return` (pbr style)
     * is used verbatim. */
    const char *sp = src;
    while (*sp == ' ' || *sp == '\t') sp++;
    int has_ret = strncmp(sp, "return", 6) == 0
                  && !(isalnum((unsigned char)sp[6]) || sp[6] == '_');
    char *wrapped = (char *)malloc(strlen(src) + 16);
    if (wrapped == NULL) { free(t); return LUA_NOREF; }
    if (has_ret) strcpy(wrapped, src);
    else sprintf(wrapped, "return (%s)", src);
    if (luaL_loadstring(L, wrapped) != LUA_OK) {
        fprintf(stderr, "[countmut] invalid %s expression: %s\n", what,
                lua_tostring(L, -1));
        lua_pop(L, 1);
        free(wrapped);
        free(t);
        return LUA_NOREF;
    }
    free(wrapped);
    free(t);
    return luaL_ref(L, LUA_REGISTRYINDEX);
}

cm_expr *cm_expr_new(const char *read_expr, const char *pile_expr) {
    if ((read_expr == NULL || read_expr[0] == '\0')
        && (pile_expr == NULL || pile_expr[0] == '\0'))
        return NULL;
    cm_expr *x = (cm_expr *)calloc(1, sizeof(*x));
    if (x == NULL) return NULL;
    lua_State *L = luaL_newstate();
    if (L == NULL) { free(x); return NULL; }
    luaL_openlibs(L);
    x->L = L;
    x->read_ref = LUA_NOREF;
    x->pile_ref = LUA_NOREF;

    lua_newtable(L); lua_pushvalue(L, -1); lua_setglobal(L, "read");
    x->read_reg = luaL_ref(L, LUA_REGISTRYINDEX);
    lua_newtable(L); lua_pushvalue(L, -1); lua_setglobal(L, "pile");
    x->pile_reg = luaL_ref(L, LUA_REGISTRYINDEX);

    /* symbolic SAM flag constants */
    struct { const char *n; int v; } flags[] = {
        {"PAIRED",1},{"PROPER_PAIR",2},{"UNMAP",4},{"MUNMAP",8},{"REVERSE",16},
        {"MREVERSE",32},{"READ1",64},{"READ2",128},{"SECONDARY",256},
        {"QCFAIL",512},{"DUP",1024},{"SUPPLEMENTARY",2048},
    };
    for (size_t i = 0; i < sizeof(flags) / sizeof(flags[0]); ++i) {
        lua_pushinteger(L, flags[i].v);
        lua_setglobal(L, flags[i].n);
    }

    lua_pushcfunction(L, l_tag);    lua_setglobal(L, "tag");
    lua_pushcfunction(L, l_exists); lua_setglobal(L, "exists");
    /* samtools-style string / regex helpers + pbr N-proportion helpers */
    lua_pushcfunction(L, l_re_match); lua_setglobal(L, "re_match");
    lua_pushcfunction(L, l_re_match); lua_setglobal(L, "re_find");
    lua_pushcfunction(L, l_slen); lua_setglobal(L, "slen");
    lua_pushcfunction(L, l_smin); lua_setglobal(L, "smin");
    lua_pushcfunction(L, l_smax); lua_setglobal(L, "smax");
    lua_pushcfunction(L, l_savg); lua_setglobal(L, "savg");
    /* htslib-style aliases: length/min/max/avg(STRING) = raw string stats */
    lua_pushcfunction(L, l_slen); lua_setglobal(L, "length");
    lua_pushcfunction(L, l_smin); lua_setglobal(L, "min");
    lua_pushcfunction(L, l_smax); lua_setglobal(L, "max");
    lua_pushcfunction(L, l_savg); lua_setglobal(L, "avg");
    /* htslib-style math helpers */
    lua_pushcfunction(L, l_sqrt); lua_setglobal(L, "sqrt");
    lua_pushcfunction(L, l_log);  lua_setglobal(L, "log");
    lua_pushcfunction(L, l_exp);  lua_setglobal(L, "exp");
    lua_pushcfunction(L, l_pow);  lua_setglobal(L, "pow");
    lua_pushcfunction(L, l_n5);   lua_setglobal(L, "n5");
    lua_pushcfunction(L, l_n3);   lua_setglobal(L, "n3");
    lua_getglobal(L, "read");
    lua_pushcfunction(L, l_n5); lua_setfield(L, -2, "n_proportion_5_prime");
    lua_pushcfunction(L, l_n3); lua_setfield(L, -2, "n_proportion_3_prime");
    lua_pop(L, 1);
    lua_pushnil(L);                 /* ensure no stale "read"/"pile" from a chunk */
    lua_pop(L, 1);

    x->read_ref = compile_chunk(
        L, read_expr, "read",
        &x->need_seq, &x->need_qual, &x->need_library,
        &x->need_qname, &x->need_rname, &x->need_mrname,
        &x->need_mapq, &x->need_flag, &x->need_strand,
        &x->need_pos, &x->need_endpos, &x->need_tid,
        &x->need_mate, &x->need_tlen, &x->need_qlen,
        &x->need_rlen, &x->need_cigar, &x->need_flags2,
        &x->need_qpos, &x->need_bq, &x->need_base,
        &x->need_ref, &x->need_dist, &x->need_aux,
        &x->need_read_tbl, &x->need_pile_tbl);
    x->pile_ref = compile_chunk(
        L, pile_expr, "pile",
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    return x;
}

int cm_expr_valid(const char *read_expr, const char *pile_expr) {
    cm_expr *x = cm_expr_new(read_expr, pile_expr);
    if (x == NULL) return 1;
    int ok
        = ((read_expr == NULL || read_expr[0] == '\0' || x->read_ref != LUA_NOREF)
           && (pile_expr == NULL || pile_expr[0] == '\0' || x->pile_ref != LUA_NOREF));
    cm_expr_free(x);
    return ok;
}

void cm_expr_free(cm_expr *x) {
    if (x == NULL) return;
    if (x->L) lua_close(x->L);
    free(x);
}

int cm_expr_has_read(const cm_expr *x) { return x != NULL && x->read_ref != LUA_NOREF; }
int cm_expr_has_pile(const cm_expr *x) { return x != NULL && x->pile_ref != LUA_NOREF; }

int cm_expr_read_constant(const cm_expr *x) {
    return x != NULL && x->read_ref != LUA_NOREF
        && !x->need_qpos && !x->need_bq && !x->need_base && !x->need_ref && !x->need_dist;
}

int cm_expr_read(cm_expr *x, const bam1_t *b, const char *rname, const char *mrname,
                 int qpos, int strand_sign, char ref_base) {
    if (x == NULL || x->read_ref == LUA_NOREF) return 1;
    lua_State *L = x->L;

    if (x->need_aux) {  /* only tag()/exists()/n5()/n3() read the current record */
        lua_pushlightuserdata(L, (void *)b);
        lua_setfield(L, LUA_REGISTRYINDEX, CUR_READ_KEY);
    }

    int lq = (int)b->core.l_qseq;
    int q = (qpos >= 0 && qpos < lq && b->core.l_qseq > 0) ? (int)bam_get_qual(b)[qpos] : -1;
    int rev = is_rev(b);
    if (qpos < 0) qpos = 0;
    int dist5 = rev ? (lq - qpos) : qpos;
    int dist3 = rev ? qpos : (lq - qpos);
    /* These are the expensive per-eval computations; run each ONLY when the
     * expression actually references a field that needs it (cigar walks and
     * the <=l_qseq base decode dominate the per-eval cost otherwise). */
    int rlen = 0;
    if (x->need_rlen || x->need_endpos) ref_span(b, &rlen);
    int sclen = 0, hclen = 0, n_indel = 0, soft5 = 0, soft3 = 0;
    if (x->need_cigar) cigar_stats(b, &sclen, &hclen, &n_indel, &soft5, &soft3);
    char seqbuf[1024], basec[2] = {'?', 0};
    const char *seq = "";
    if (x->need_seq || x->need_base) {
        seq = read_seq_str(b, seqbuf, sizeof(seqbuf));
        if (x->need_base && qpos >= 0 && (size_t)qpos < (size_t)lq) basec[0] = seq[qpos];
    }
    char refc[2] = { ref_base ? ref_base : 'N', 0 };

    /* mate position (1-based, 0 if none); insert size (TLEN) */
    int mpos = (b->core.mpos >= 0) ? b->core.mpos + 1 : 0;
    int tlen = (int)b->core.isize;

    if (x->need_qname)   r_str(x, "qname", bam_get_qname(b));
    if (x->need_rname)   r_str(x, "rname", rname ? rname : "");
    if (x->need_mrname) {
        r_str(x, "mrname", mrname ? mrname : "");
        r_str(x, "rnext", mrname ? mrname : "");
    }
    if (x->need_seq) {
        r_str(x, "seq", seq);
        r_str(x, "sequence", seq);
    }
    if (x->need_qual) {
        char *qb = seqbuf;           /* reuse the buffer */
        size_t qn = (size_t)lq;
        if (qn > 1023) qn = 1023;
        for (size_t z = 0; z < qn; ++z) qb[z] = (char)bam_get_qual(b)[z];
        qb[qn] = 0;
        r_str(x, "qual", qb);
    }
    if (x->need_library) {
        uint8_t *lb = bam_aux_get(b, "LB");
        r_str(x, "library", (lb && *lb == 'Z') ? (const char *)(lb + 1) : "");
    }
    if (x->need_base) r_str(x, "base", basec);
    if (x->need_ref)  r_str(x, "ref", refc);
    int strandv = strand_sign;        /* +1 forward / -1 reverse (biological) */
    if (x->need_mapq)   { r_int(x, "mapq", (int)b->core.qual); r_int(x, "MAPQ", (int)b->core.qual); r_int(x, "mapping_quality", (int)b->core.qual); }
    if (x->need_flag)   { r_int(x, "flag", (int)b->core.flag); r_int(x, "flags", (int)b->core.flag); r_int(x, "FLAGS", (int)b->core.flag); }
    if (x->need_strand) { r_int(x, "strand", strandv); r_int(x, "STRAND", strandv); }
    if (x->need_pos)    { r_int(x, "pos", (int)b->core.pos + 1); r_int(x, "POS", (int)b->core.pos + 1); r_int(x, "start", (int)b->core.pos + 1); }
    if (x->need_endpos) { r_int(x, "endpos", (int)b->core.pos + rlen + 1); r_int(x, "stop", (int)b->core.pos + rlen + 1); }
    if (x->need_tid)    { r_int(x, "tid", b->core.tid); r_int(x, "refid", b->core.tid); }
    if (x->need_mate)   { r_int(x, "mtid", b->core.mtid); r_int(x, "mpos", mpos); r_int(x, "pnext", mpos); }
    if (x->need_tlen)   { r_int(x, "tlen", tlen); r_int(x, "insert_size", tlen); }
    if (x->need_qpos) {
        r_int(x, "qpos", qpos);
        r_int(x, "QPOS", qpos);
    }
    if (x->need_qlen)   { r_int(x, "qlen", lq); r_int(x, "length", lq); r_int(x, "LEN", lq); }
    if (x->need_rlen)   r_int(x, "rlen", rlen);
    if (x->need_cigar) {
        r_int(x, "sclen", sclen);
        r_int(x, "hclen", hclen);
        r_int(x, "n_indel", n_indel);
        r_int(x, "indel_count", n_indel);
        r_int(x, "soft_clips_5_prime", soft5);
        r_int(x, "soft_clips_3_prime", soft3);
        r_int(x, "ncigar", (int)b->core.n_cigar);
    }
    if (x->need_bq) {
        r_int(x, "bq", q);
        r_int(x, "BQ", q);
        r_int(x, "baseq", q);
    }
    if (x->need_flags2) {
        r_int(x, "is_reverse", rev ? 1 : 0);
        r_int(x, "is_paired", (b->core.flag & BAM_FPAIRED) ? 1 : 0);
        r_int(x, "r1", (b->core.flag & BAM_FREAD1) ? 1 : 0);
        r_int(x, "r2", (b->core.flag & 128) ? 1 : 0);
    }
    if (x->need_dist) {
        r_int(x, "distance_from_5prime", dist5);
        r_int(x, "dist5", dist5);
        r_int(x, "DIST5", dist5);
        r_int(x, "distance_from_3prime", dist3);
        r_int(x, "dist3", dist3);
        r_int(x, "DIST3", dist3);
    }

    return run_chunk(L, x->read_ref);
}

int cm_expr_pile(cm_expr *x, int64_t pos, char ref_ch, const char *motif,
                 const int cnt[5], int ins, int del, int rs, int fl) {
    if (x == NULL || x->pile_ref == LUA_NOREF) return 1;
    lua_State *L = x->L;
    int depth = 0;
    for (int i = 0; i < 5; ++i) depth += cnt[i];
    char refstr[2] = { ref_ch ? ref_ch : 'N', 0 };
    const char *mot = motif ? motif : "";

    p_int(x, "pos", (lua_Integer)pos);
    p_int(x, "POS", (lua_Integer)pos);
    p_str(x, "ref", refstr);
    p_str(x, "ref_base", refstr);
    p_str(x, "REF", refstr);
    p_int(x, "depth", depth);
    p_int(x, "DEPTH", depth);
    p_int(x, "a", cnt[0]);
    p_int(x, "c", cnt[1]);
    p_int(x, "g", cnt[2]);
    p_int(x, "t", cnt[3]);
    p_int(x, "n", cnt[4]);
    p_int(x, "A", cnt[0]);
    p_int(x, "C", cnt[1]);
    p_int(x, "G", cnt[2]);
    p_int(x, "T", cnt[3]);
    p_int(x, "N", cnt[4]);
    p_int(x, "ins", ins);
    p_int(x, "del", del);
    p_int(x, "ref_skip", rs);
    p_int(x, "fail", fl);
    p_str(x, "motif", mot);

    return run_chunk(L, x->pile_ref);
}
