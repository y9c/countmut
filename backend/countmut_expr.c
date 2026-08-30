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

#include "countmut_expr.h"
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#define CUR_READ_KEY "countmut.request_read"

struct cm_expr {
    lua_State *L;
    int read_ref;   /* LUA_NOREF if no -e */
    int pile_ref;   /* LUA_NOREF if no -p */
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

static int is_rev(const bam1_t *b) { return (b->core.flag & 16) != 0; }

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

/* set one scalar into the global namespace AND into the `read`/`pile` table */
static void put_int(lua_State *L, const char *tbl, const char *name, lua_Integer v) {
    lua_pushinteger(L, v); lua_setglobal(L, name);
    lua_getglobal(L, tbl);
    lua_pushinteger(L, v); lua_setfield(L, -2, name);
    lua_pop(L, 1);
}
static void put_str(lua_State *L, const char *tbl, const char *name, const char *v) {
    lua_pushstring(L, v); lua_setglobal(L, name);
    lua_getglobal(L, tbl);
    lua_pushstring(L, v); lua_setfield(L, -2, name);
    lua_pop(L, 1);
}

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

/* countmut/Python/samtools-style expressions use !=, &&, || and !; Lua uses
 * ~=, and, or, not.  Translate (outside string literals) so existing -e / -p
 * expressions keep working. */
static char *translate_ops(const char *src) {
    size_t n = strlen(src);
    char *out = (char *)malloc(2 * n + 1);
    if (out == NULL) return NULL;
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
        out[j++] = c;
    }
    out[j] = 0;
    return out;
}

static int compile_chunk(lua_State *L, const char *expr, const char *what) {
    if (expr == NULL || expr[0] == '\0') return LUA_NOREF;
    char *t = translate_ops(expr);
    const char *src = t ? t : expr;
    if (luaL_loadstring(L, src) != LUA_OK) {
        /* Lua only accepts function-call statements as bare expression chunks,
         * so pbr-style has to write `return ...`.  countmut's documented -e/-p
         * are *bare predicates* (`mapq >= 20`); wrap them in return (...). */
        lua_pop(L, 1);
        char *wrapped = (char *)malloc(strlen(src) + 16);
        if (wrapped == NULL) { free(t); return LUA_NOREF; }
        sprintf(wrapped, "return (%s)", src);
        int ok = (luaL_loadstring(L, wrapped) == LUA_OK);
        if (!ok) {
            fprintf(stderr, "[countmut] invalid %s expression: %s\n", what,
                    lua_tostring(L, -1));
            lua_pop(L, 1);
        }
        free(wrapped);
        free(t);
        if (!ok) return LUA_NOREF;
        return luaL_ref(L, LUA_REGISTRYINDEX);
    }
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

    lua_newtable(L); lua_setglobal(L, "read");
    lua_newtable(L); lua_setglobal(L, "pile");

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
    lua_pushnil(L);                 /* ensure no stale "read"/"pile" from a chunk */
    lua_pop(L, 1);

    x->read_ref = compile_chunk(L, read_expr, "read");
    x->pile_ref = compile_chunk(L, pile_expr, "pile");
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

int cm_expr_read(cm_expr *x, const bam1_t *b, const char *rname, int qpos, int strand_sign) {
    if (x == NULL || x->read_ref == LUA_NOREF) return 1;
    lua_State *L = x->L;

    lua_pushlightuserdata(L, (void *)b);
    lua_setfield(L, LUA_REGISTRYINDEX, CUR_READ_KEY);

    int lq = (int)b->core.l_qseq;
    int q = (qpos >= 0 && qpos < lq && b->core.l_qseq > 0) ? (int)bam_get_qual(b)[qpos] : -1;
    int rev = is_rev(b);
    if (qpos < 0) qpos = 0;
    int dist5 = rev ? (lq - qpos) : qpos;
    int dist3 = rev ? qpos : (lq - qpos);
    int rlen = 0; ref_span(b, &rlen);

    put_str(L, "read", "qname", bam_get_qname(b));
    put_str(L, "read", "rname", rname ? rname : "");
    int strandv = strand_sign;        /* +1 forward / -1 reverse (biological) */
    put_int(L, "read", "mapq", (int)b->core.qual);
    put_int(L, "read", "MAPQ", (int)b->core.qual);
    put_int(L, "read", "flag", (int)b->core.flag);
    put_int(L, "read", "flags", (int)b->core.flag);
    put_int(L, "read", "FLAGS", (int)b->core.flag);
    put_int(L, "read", "strand", strandv);
    put_int(L, "read", "STRAND", strandv);
    put_int(L, "read", "pos", (int)b->core.pos + 1);
    put_int(L, "read", "POS", (int)b->core.pos + 1);
    put_int(L, "read", "endpos", (int)b->core.pos + rlen + 1);
    put_int(L, "read", "qpos", qpos);
    put_int(L, "read", "QPOS", qpos);
    put_int(L, "read", "qlen", lq);
    put_int(L, "read", "length", lq);
    put_int(L, "read", "LEN", lq);
    put_int(L, "read", "rlen", rlen);
    put_int(L, "read", "bq", q);
    put_int(L, "read", "BQ", q);
    put_int(L, "read", "baseq", q);
    put_int(L, "read", "distance_from_5prime", dist5);
    put_int(L, "read", "dist5", dist5);
    put_int(L, "read", "DIST5", dist5);
    put_int(L, "read", "distance_from_3prime", dist3);
    put_int(L, "read", "dist3", dist3);
    put_int(L, "read", "DIST3", dist3);

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

    put_int(L, "pile", "pos", (lua_Integer)pos);
    put_int(L, "pile", "POS", (lua_Integer)pos);
    put_str(L, "pile", "ref", refstr);
    put_str(L, "pile", "ref_base", refstr);
    put_str(L, "pile", "REF", refstr);
    put_int(L, "pile", "depth", depth);
    put_int(L, "pile", "DEPTH", depth);
    put_int(L, "pile", "a", cnt[0]);
    put_int(L, "pile", "c", cnt[1]);
    put_int(L, "pile", "g", cnt[2]);
    put_int(L, "pile", "t", cnt[3]);
    put_int(L, "pile", "n", cnt[4]);
    put_int(L, "pile", "A", cnt[0]);
    put_int(L, "pile", "C", cnt[1]);
    put_int(L, "pile", "G", cnt[2]);
    put_int(L, "pile", "T", cnt[3]);
    put_int(L, "pile", "N", cnt[4]);
    put_int(L, "pile", "ins", ins);
    put_int(L, "pile", "del", del);
    put_int(L, "pile", "ref_skip", rs);
    put_int(L, "pile", "fail", fl);
    put_str(L, "pile", "motif", mot);

    return run_chunk(L, x->pile_ref);
}
