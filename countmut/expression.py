#!/usr/bin/env python3
"""
Samtools-style filter-expression engine for countmut.

Implements the *samtools ``--input-fmt-option filter=STRING``* expression
grammar (C-style precedence, ``&&``/``||``/``!``, bit fields ``flag.dup``,
tag refs ``[NM]``, regex ``=~``/``!~``, and the samtools SAM-field variable
namespace) plus a small site-level extension for ``-p``.

Two string filters are provided:

* ``-e, --expression <STR>`` -- per-aligned-base *read* filter (samtools vars).
* ``-p, --pile-expression <STR>`` -- per-site filter (pileup vars).

Author: Ye Chang
Date: 2026-08-30
"""

from __future__ import annotations

import math
import re

import pysam

from . import reads

# ---------------------------------------------------------------------------
# SAM flag bits (samtools naming) for flag.<name>
# ---------------------------------------------------------------------------
FLAG_BITS = {
    "paired": 1,
    "proper_pair": 2,
    "unmap": 4,
    "munmap": 8,
    "reverse": 16,
    "mreverse": 32,
    "read1": 64,
    "read2": 128,
    "secondary": 256,
    "qcfail": 512,
    "dup": 1024,
    "supplementary": 2048,
}
FLAG_BITS["proper-pair"] = 2
FLAG_BITS["read-1"] = 64
FLAG_BITS["read-2"] = 128

# ---------------------------------------------------------------------------
# Lexer
# ---------------------------------------------------------------------------
_OP3 = ("=~", "!~", "...")
_OP2 = ("&&", "||", "==", "!=", "<=", ">=")
_OP1 = set("!~&^|+-*/%<>()[],.")


def _lex(src: str):
    toks = []
    i, n = 0, len(src)
    while i < n:
        c = src[i]
        if c.isspace():
            i += 1
            continue
        # number (decimal / float / 0x hex)
        if c.isdigit() or (c == "." and i + 1 < n and src[i + 1].isdigit()):
            j = i
            if src.startswith("0x", i) or src.startswith("0X", i):
                j = i + 2
                while j < n and (src[j].isdigit() or src[j].lower() in "abcdef"):
                    j += 1
                toks.append(("num", int(src[i + 2 : j], 16)))
                i = j
                continue
            while j < n and src[j].isdigit():
                j += 1
            if j < n and src[j] == ".":
                j += 1
                while j < n and src[j].isdigit():
                    j += 1
                toks.append(("num", float(src[i:j])))
                i = j
                continue
            toks.append(("num", int(src[i:j])))
            i = j
            continue
        # string
        if c in ("'", '"'):
            q = c
            i += 1
            buf = []
            while i < n and src[i] != q:
                if src[i] == "\\" and i + 1 < n:
                    buf.append(src[i + 1])
                    i += 2
                else:
                    buf.append(src[i])
                    i += 1
            i += 1  # closing quote
            toks.append(("str", "".join(buf)))
            continue
        # tag ref  [XX]
        if c == "[":
            j = src.find("]", i)
            if j == -1:
                raise ValueError("unterminated tag reference")
            toks.append(("tag", src[i + 1 : j]))
            i = j + 1
            continue
        # identifier (also accept and/or/not as aliases for &&/||/!)
        if c.isalpha() or c == "_":
            j = i
            while j < n and (src[j].isalnum() or src[j] in "_:-"):
                j += 1
            word = src[i:j]
            if word == "and":
                toks.append(("op", "&&"))
            elif word == "or":
                toks.append(("op", "||"))
            elif word == "not" or word == "!":
                toks.append(("op", "!"))
            else:
                toks.append(("ident", word))
            i = j
            continue
        # operators (longest match)
        three = src[i : i + 3]
        two = src[i : i + 2]
        if three in _OP3:
            toks.append(("op", three))
            i += 3
            continue
        if two in _OP2:
            toks.append(("op", two))
            i += 2
            continue
        if c in _OP1:
            toks.append(("op", c))
            i += 1
            continue
        raise ValueError(f"unexpected character '{c}' in expression")
    toks.append(("eof", None))
    return toks


# ---------------------------------------------------------------------------
# Parser (precedence climbing, samtools precedence)
# ---------------------------------------------------------------------------
# binding power of a binary operator (centre, left)
_BINFIX = {
    "||": (70, 71),
    "&&": (80, 81),
    "==": (90, 91),
    "!=": (90, 91),
    "=~": (90, 91),
    "!~": (90, 91),
    ">": (100, 101),
    ">=": (100, 101),
    "<": (100, 101),
    "<=": (100, 101),
    "|": (110, 111),
    "^": (120, 121),
    "&": (130, 131),
    "+": (140, 141),
    "-": (140, 141),
    "*": (150, 151),
    "/": (150, 151),
    "%": (150, 151),
}
_UNARY_BP = 160


class _Parser:
    def __init__(self, toks):
        self.toks = toks
        self.i = 0

    def peek(self):
        return self.toks[self.i]

    def next(self):
        t = self.toks[self.i]
        self.i += 1
        return t

    def parse(self):
        node = self.expr(0)
        if self.peek()[0] != "eof":
            raise ValueError("trailing tokens in expression")
        return node

    def expr(self, min_bp):
        # prefix
        tok = self.next()
        kind, val = tok
        if kind == "op" and val in ("!", "~", "-", "+"):
            right = self.expr(_UNARY_BP)
            node = ("unary", val, right)
        else:
            node = self.primary(tok)
        # postfix function calls (bind tightest):  name(a, b)
        while self.peek() == ("op", "("):
            self.next()
            args = []
            if self.peek() != ("op", ")"):
                args.append(self.expr(0))
                while self.peek() == ("op", ","):
                    self.next()
                    args.append(self.expr(0))
            if self.peek() == ("op", ")"):
                self.next()
            node = ("call", node, args)
        # infix
        while True:
            kind, val = self.peek()
            if kind != "op" or val not in _BINFIX:
                break
            lbp, rbp = _BINFIX[val]
            if lbp < min_bp:
                break
            self.next()
            right = self.expr(rbp)
            node = ("binop", val, node, right)
        return node

    def primary(self, tok):
        kind, val = tok
        if kind in ("num", "str"):
            return ("lit", val)
        if kind == "ident":
            # flag.<bit> dotted fields
            if val.lower() == "flag" and self.peek() == ("op", "."):
                self.next()
                bit = self.next()
                if bit[0] == "ident":
                    return ("flagbit", bit[1])
                raise ValueError("expected flag bit name after 'flag.'")
            if self.peek() == ("op", "."):
                # generic dotted access -> attribute of a value (rare)
                self.next()
                attr = self.next()
                if attr[0] == "ident":
                    return ("attr", val, attr[1])
            return ("var", val)
        if kind == "tag":
            return ("tag", val)
        if kind == "op" and val == "(":
            node = self.expr(0)
            if self.peek() == ("op", ")"):
                self.next()
            return node
        raise ValueError(f"unexpected token {kind}:{val} in expression")


def _parse(src):
    return _Parser(_lex(src)).parse()


# ---------------------------------------------------------------------------
# Evaluator
# ---------------------------------------------------------------------------
def _truthy(v):
    return bool(v)  # None -> False


def _eval(node, ns, funcs):
    kind = node[0]
    if kind == "lit":
        return node[1]
    if kind == "var":
        return ns.get(node[1])  # None if undefined
    if kind == "tag":
        return ns.get("__tag__")(node[1]) if "__tag__" in ns else None
    if kind == "flagbit":
        flag = ns.get("flag", 0)
        return flag & FLAG_BITS.get(node[1], 0)
    if kind == "attr":
        v = _eval(node[1], ns, funcs)
        return getattr(v, node[2], None) if v is not None else None
    if kind == "call":
        fn = _eval(node[1], ns, funcs)
        args = [_eval(a, ns, funcs) for a in node[2]]
        return fn(*args) if callable(fn) else None
    if kind == "unary":
        v = _eval(node[2], ns, funcs)
        op = node[1]
        if op == "!":
            return not _truthy(v)
        if op == "-":
            return -v
        if op == "+":
            return +v
        if op == "~":
            return ~int(v)
    if kind == "binop":
        op, a, b = node[1], _eval(node[2], ns, funcs), _eval(node[3], ns, funcs)
        return _binop(op, a, b, ns, funcs)
    raise ValueError(f"unknown node {kind}")


def _binop(op, a, b, ns, funcs):
    if op == "&&":
        return _truthy(a) and _truthy(b)
    if op == "||":
        return _truthy(a) or _truthy(b)
    # comparisons: undefined (None) -> false, except != / !~
    if op in ("==", "!=", "<", "<=", ">", ">=", "=~", "!~"):
        if a is None or b is None:
            return False
        if op == "==":
            return a == b
        if op == "!=":
            return a != b
        if op == "=~":
            return re.search(str(b), str(a)) is not None
        if op == "!~":
            return re.search(str(b), str(a)) is None
        if op == ">":
            return a > b
        if op == ">=":
            return a >= b
        if op == "<":
            return a < b
        if op == "<=":
            return a <= b
    if op == "+":
        return a + b
    if op == "-":
        return a - b
    if op == "*":
        return a * b
    if op == "/":
        return a / b
    if op == "%":
        return int(a) % int(b) if b else 0
    if op == "&":
        return int(a) & int(b)
    if op == "^":
        return int(a) ^ int(b)
    if op == "|":
        return int(a) | int(b)
    return False


def _apply_func(name, args):
    if args is None or len(args) == 0:
        return None
    if name in ("length", "len"):
        return len(args[0]) if args[0] is not None else 0
    if name == "min":
        return min(_as_iter(args[0])) if args[0] is not None else None
    if name == "max":
        return max(_as_iter(args[0])) if args[0] is not None else None
    if name == "avg":
        it = _as_iter(args[0]) if args[0] is not None else []
        return sum(it) / len(it) if it else float("nan")
    if name == "exists":
        return _truthy(args[0])
    if name == "default":
        return args[0] if _truthy(args[0]) else (args[1] if len(args) > 1 else None)
    if name == "sqrt":
        return math.sqrt(args[0])
    if name == "log":
        return math.log(args[0])
    if name == "exp":
        return math.exp(args[0])
    if name == "pow":
        return math.pow(args[0], args[1])
    return None


def _fn(name):
    """Bind a function NAME so that ``name(x)`` dispatches correctly."""
    return lambda *a: _apply_func(name, a)


def _as_iter(v):
    if isinstance(v, (bytes, str)):
        return list(v)
    if v is None:
        return []
    return list(v)


# ---------------------------------------------------------------------------
# public compile helpers
# ---------------------------------------------------------------------------
def _compile(src, stack):
    ast = _parse(src)
    expr = ast  # node

    def run(ns, funcs):
        try:
            return bool(_eval(expr, ns, funcs))
        except Exception:
            return False

    return run


COMPILED = {}


def compile_read_pred(src):
    """samtools SAM-field read predicate -> ``pred(read, qpos)``."""
    key = ("r", src)
    if key not in COMPILED:
        COMPILED[key] = _compile(src, None)
    run = COMPILED[key]

    def pred(rec: pysam.AlignedSegment, qpos):
        ns = _read_ns(rec, qpos)
        funcs = {"tag": rec.get_tag}
        return run(ns, funcs)

    return pred


def compile_pile_pred(src):
    """pileup-site predicate -> ``pred(site_column)``."""
    key = ("p", src)
    if key not in COMPILED:
        COMPILED[key] = _compile(src, None)
    run = COMPILED[key]

    def pred(col):
        ns = _pile_ns(col)
        return run(ns, {})

    return pred


def _read_ns(rec, qpos):
    cig = rec.cigartuples or []
    rlen = 0
    for op, ln in cig:
        if op in (0, 2, 3, 7, 8):
            rlen += ln
    sclen = sum(ln for op, ln in cig if op == 4)
    hclen = sum(ln for op, ln in cig if op == 5)
    q = rec.query_qualities
    bq = int(q[qpos]) if (q is not None and qpos is not None and qpos < len(q)) else -1
    length = len(rec.query_sequence or "")
    strand = 1 if reads.actual_strand(rec) == "+" else -1
    qname = rec.query_name

    def tag(name):
        try:
            return rec.get_tag(name)
        except (KeyError, ValueError):
            return None

    return {
        # symbolic flag constants (samtools style): flag & UNMAP == 0
        **{k.upper(): v for k, v in FLAG_BITS.items()},
        "mapq": rec.mapping_quality,
        "flag": rec.flag,
        "qname": qname,
        "pos": rec.reference_start + 1,
        "endpos": rec.reference_end or (rec.reference_start + 1),
        "pnext": (rec.next_reference_start + 1) if rec.next_reference_start >= 0 else 0,
        "mpos": (rec.next_reference_start + 1) if rec.next_reference_start >= 0 else 0,
        "rname": rec.reference_name or "",
        "mrname": rec.next_reference_name or "",
        "tlen": rec.template_length,
        "qlen": length,
        "rlen": rlen,
        "ncigar": len(cig),
        "seq": rec.query_sequence or "",
        "qual": q if q is not None else "",
        "sclen": sclen,
        "hclen": hclen,
        "strand": strand,
        "bq": bq,
        "distance_from_5prime": (length - qpos)
        if (qpos is not None and rec.is_reverse)
        else (qpos or 0),
        "distance_from_3prime": qpos
        if (qpos is not None and rec.is_reverse)
        else (length - (qpos or 0)),
        "library": tag("LB"),
        "tag": tag,
        "__tag__": tag,
        "exists": _apply_func_exists,
        "default": _apply_func_default,
        "length": _fn("length"),
        "min": _fn("min"),
        "max": _fn("max"),
        "avg": _fn("avg"),
        "sqrt": _fn("sqrt"),
        "log": _fn("log"),
        "pow": _fn("pow"),
        "exp": _fn("exp"),
    }


def _pile_ns(col):
    a = sum(col.counts.get(s, {}).get("base", {}).get("A", 0) for s in ("+", "-"))
    c = sum(col.counts.get(s, {}).get("base", {}).get("C", 0) for s in ("+", "-"))
    g = sum(col.counts.get(s, {}).get("base", {}).get("G", 0) for s in ("+", "-"))
    t = sum(col.counts.get(s, {}).get("base", {}).get("T", 0) for s in ("+", "-"))
    n = sum(col.counts.get(s, {}).get("base", {}).get("N", 0) for s in ("+", "-"))
    return {
        "depth": col.total_depth(),
        "pos": col.pos,
        "ref": col.ref_base,
        "ref_base": col.ref_base,
        "a": a,
        "c": c,
        "g": g,
        "t": t,
        "n": n,
        "A": a,
        "C": c,
        "G": g,
        "T": t,
        "N": n,
        "ins": col.ins.get("+", 0) + col.ins.get("-", 0),
        "del": col.deletes.get("+", 0) + col.deletes.get("-", 0),
        "ref_skip": col.ref_skip.get("+", 0) + col.ref_skip.get("-", 0),
        "fail": col.fail.get("+", 0) + col.fail.get("-", 0),
        "motif": col.motif,
        "exists": _apply_func_exists,
        "default": _apply_func_default,
        "length": _fn("length"),
        "min": _fn("min"),
        "max": _fn("max"),
        "avg": _fn("avg"),
        "sqrt": _fn("sqrt"),
        "log": _fn("log"),
        "pow": _fn("pow"),
        "exp": _fn("exp"),
    }


def _apply_func_exists(*args):
    return _truthy(args[0]) if args else False


def _apply_func_default(*args):
    return (
        args[0]
        if (_truthy(args[0]) if args else False)
        else (args[1] if len(args) > 1 else None)
    )
