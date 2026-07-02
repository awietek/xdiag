# SPDX-License-Identifier: Apache-2.0
#
# Emit the ergonomic Julia side of the wrapper from the API model + overrides.
#
# For every wrapped C++ type T (registered on the C++ side as cxx_T) this emits
# a Julia struct T holding a cxx_T, a convert() from cxx_T, constructor
# forwarders (with Julia default arguments and 1-based index shifts), and method
# forwarders that unwrap the cxx_ fields of their arguments, call the CxxWrap
# function, and re-wrap / shift the result. Blocks are handled through the
# abstract Julia type Block (the concrete cxx block dispatches on the C++ side).
#
# The C++ side (emit_cpp.py) is a thin pass-through; all the 1-based<->0-based
# shifts and default arguments live here.

import itertools
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import model as _model
import overrides as ov

BLOCK_TYPES = {"Spinhalf", "tJ", "Electron", "Boson", "Fermion"}
RESULT_STRUCTS = {
    "EigsLanczosResult", "EigvalsLanczosResult", "EvolveLanczosResult",
    "EvolveLanczosInplaceResult", "TimeEvolveExpokitResult",
    "TimeEvolveExpokitInplaceResult", "EigsLobpcgResult",
}
ITERATORS = {  # iterator type -> element is a ProductState
    "SpinhalfIterator", "tJIterator", "ElectronIterator", "BosonIterator",
    "FermionIterator",
}

# Method names that extend Base / LinearAlgebra (must be imported in XDiag.jl).
BASE_NAMES = {
    "+", "-", "*", "/", "==", "!=", "getindex", "setindex!", "size", "isreal",
    "show", "real", "imag", "push!", "length", "isapprox", "inv", "adjoint",
    "conj", "abs",
}
LINALG_NAMES = {"norm", "dot"}

# C++ default value -> Julia default value.
_DEFAULT_MAP = {"true": "true", "false": "false", "nullptr": "nothing"}


def short(q):
    return q.split("::")[-1]


def _snake(name):
    out = []
    for i, c in enumerate(name):
        if c.isupper() and i > 0 and not name[i - 1].isupper():
            out.append("_")
        out.append(c.lower())
    return "".join(out)


class Skip(Exception):
    pass


def field_name(qual):
    return ov.FIELD_NAMES.get(qual, "cxx_" + _snake(short(qual)))


# ---------------------------------------------------------------------------
# Type mapping: Julia parameter type, how to unwrap an arg, how to wrap a return
# ---------------------------------------------------------------------------

def jl_param_type(tp, wrapped):
    cat = tp["cat"]
    if cat == "primitive":
        return {"bool": "Bool", "int64": "Int64", "int32": "Int32",
                "double": "Float64", "complex": "ComplexF64",
                "uint64": "UInt64", "uint32": "UInt32"}.get(tp.get("prim"), "")
    if cat == "string":
        return "String"
    if cat == "variant" and tp.get("is_block"):
        return "Block"
    if cat == "wrapped":
        q = tp["type"]
        if q in ov.NATIVE_BRIDGES:
            return {"scalar": "Number", "coeff": "Union{Number,String}",
                    "vector": "Vector", "matrix": "Matrix"}[ov.NATIVE_BRIDGES[q]]
        if q in wrapped:
            return short(q)
        raise Skip("param wrapped " + q)
    if cat == "arma":
        return {"vec": "Vector{Float64}", "cx_vec": "Vector{ComplexF64}",
                "mat": "Matrix{Float64}", "cx_mat": "Matrix{ComplexF64}",
                "ivec": "Vector{Int64}"}.get(tp.get("arma"), "")
    if cat == "stdvector":
        e = tp.get("elem", {})
        if e.get("cat") == "primitive":
            jt = {"int64": "Int64", "int32": "Int32", "double": "Float64",
                  "complex": "ComplexF64"}.get(e.get("prim"), "")
            return f"Vector{{{jt}}}" if jt else "Vector"
        if e.get("cat") == "string":
            return "Vector{String}"
        if e.get("cat") == "wrapped":
            return f"Vector{{{short(e['type'])}}}"
        raise Skip("param vector")
    raise Skip("param " + cat)


def unwrap_arg(tp, name, wrapped, shift):
    """Julia expression passing `name` to the cxx function."""
    cat = tp["cat"]
    if cat == "variant" and tp.get("is_block"):
        return f"{name}.cxx_block"
    if cat == "wrapped":
        q = tp["type"]
        if q in ov.NATIVE_BRIDGES:
            # scalar/coeff pass through as native Julia numbers; vector/matrix
            # are marshalled through armadillo (C++ overloads take arma types).
            if ov.NATIVE_BRIDGES[q] in ("vector", "matrix"):
                return f"to_armadillo({name})"
            return name
        return f"{name}.{field_name(q)}"
    if cat == "stdvector":
        e = tp.get("elem", {})
        inner = f"{name} .- 1" if shift else name
        if e.get("cat") == "wrapped":
            fn = field_name(e["type"])
            base = f"[x.{fn} for x in {name}]"
            return f"StdVector({base})"
        return f"StdVector({inner})"
    if cat == "arma":
        return f"to_armadillo({name} .- 1)" if shift else f"to_armadillo({name})"
    if cat == "primitive" and shift:
        return f"{name} - 1"
    return name


def wrap_return(tp, expr, wrapped, shift):
    cat = tp["cat"]
    if cat == "void":
        return expr, False
    if cat == "wrapped":
        q = tp["type"]
        if q in ov.NATIVE_BRIDGES:
            # vector/matrix come back as armadillo objects; scalar/coeff are
            # already native Julia numbers.
            if ov.NATIVE_BRIDGES[q] in ("vector", "matrix"):
                return f"to_julia({expr})", True
            return expr, True
        if q not in wrapped:
            raise Skip("return wrapped " + q)
        return f"{short(q)}({expr})", True
    if cat == "arma":
        return f"to_julia({expr})", True
    if cat == "stdvector":
        e = tp.get("elem", {})
        jt = jl_param_type(tp, wrapped)
        v = f"{jt}({expr})"
        if shift:
            v = f"({v} .+ 1)"
        return v, True
    if cat == "primitive" and shift:
        return f"({expr} + 1)", True
    if cat == "variant":
        raise Skip("return variant")
    if cat == "std_composite":
        # std::tuple returns (eig0/eigs) pass through as a Julia tuple; other
        # std composites (iterators, chrono, ...) are not exposed.
        if tp.get("cpp", "").startswith("std::tuple<"):
            return expr, True
        raise Skip("return std_composite")
    if cat in ("template_param", "pointer", "unknown", "ostream", "enum"):
        raise Skip("return " + cat)
    return expr, True


def default_jl(cpp_default):
    if cpp_default is None:
        return None
    d = cpp_default.strip()
    if d in _DEFAULT_MAP:
        return _DEFAULT_MAP[d]
    if d.startswith('"'):
        return d  # string literal
    d = d.rstrip("Ll")  # 1000L -> 1000
    return d  # numeric literal passes through


# ---------------------------------------------------------------------------
# Emission of one callable (constructor / method / free function)
# ---------------------------------------------------------------------------

def onebased_args(qual):
    return set(ov.ONEBASED_ARGS.get(qual, []))


_JL_KEYWORDS = {
    "begin", "end", "function", "if", "else", "elseif", "while", "for", "do",
    "let", "quote", "module", "baremodule", "struct", "mutable", "return",
    "in", "isa", "where", "try", "catch", "finally", "const", "global", "local",
    "import", "using", "export", "abstract", "primitive", "type", "macro",
}


def _safe_name(nm):
    return nm + "_" if nm in _JL_KEYWORDS else nm


def build_params(entry, wrapped, self_type=None, self_field=None):
    """Return (param_decls, callargs, names). Raises Skip on an unmappable arg."""
    decls, callargs, names = [], [], []
    if self_type is not None:
        decls.append(f"self::{self_type}")
    shifts = onebased_args(entry["qualified"])
    for i, a in enumerate(entry.get("args", [])):
        jt = jl_param_type(a["type"], wrapped)
        orig = a["name"] or f"a{i}"
        nm = _safe_name(orig)
        decl = f"{nm}::{jt}" if jt else nm
        dv = default_jl(a.get("default")) if a.get("has_default") else None
        if dv is not None:
            decl += f"={dv}"
        decls.append(decl)
        names.append(nm)
        callargs.append(unwrap_arg(a["type"], nm, wrapped, orig in shifts))
    return decls, callargs, names


def emit_constructor(entry, wrapped):
    rec = short(entry["record"])
    decls, callargs, _ = build_params(entry, wrapped)
    call = f"construct_{rec}({', '.join(callargs)})"
    return f"{rec}({', '.join(decls)}) = {rec}({call})"


def emit_method(entry, wrapped):
    rec = short(entry["record"])
    fld = field_name(entry["qualified"].rsplit("::", 1)[0])
    ret_shift = entry["qualified"] in ov.ONEBASED_RETURN
    decls, callargs, _ = build_params(entry, wrapped, self_type=rec)
    name = entry["name"]
    if entry["kind"] == "operator":
        raise Skip("operator (handled elsewhere)")
    call = f"{name}(self.{fld}{', ' if callargs else ''}{', '.join(callargs)})"
    wrapped_expr, _ = wrap_return(entry.get("return", {"cat": "void"}), call,
                                  wrapped, ret_shift)
    return f"{name}(" + ", ".join(decls) + f") = {wrapped_expr}"



# ---------------------------------------------------------------------------
# operators, free functions, field accessors
# ---------------------------------------------------------------------------

_JL_OP = {"+": "+", "-": "-", "*": "*", "/": "/", "==": "(==)", "!=": "(!=)"}


def emit_operator(entry, wrapped, member):
    op = entry.get("op")
    nargs = len(entry.get("args", []))
    ret = entry.get("return", {"cat": "void"})
    ret_shift = entry["qualified"] in ov.ONEBASED_RETURN
    shifts = onebased_args(entry["qualified"])
    if member:
        rec = short(entry["record"])
        fld = field_name(entry["qualified"].rsplit("::", 1)[0])
        if op in _JL_OP and nargs == 1:
            a = entry["args"][0]
            jt = jl_param_type(a["type"], wrapped)
            nm = a["name"] or "b"
            rhs = unwrap_arg(a["type"], nm, wrapped, nm in shifts)
            call = f"{op}(self.{fld}, {rhs})"
            expr, _ = wrap_return(ret, call, wrapped, ret_shift)
            return f"Base.:{_JL_OP[op]}(self::{rec}, {nm}::{jt}) = {expr}"
        if op == "[]" and nargs == 1:
            a = entry["args"][0]
            nm = a["name"] or "i"
            rhs = unwrap_arg(a["type"], nm, wrapped, nm in shifts)
            call = f"getindex(self.{fld}, {rhs})"
            expr, _ = wrap_return(ret, call, wrapped, ret_shift)
            return f"Base.getindex(self::{rec}, {nm}::Int64) = {expr}"
        raise Skip(f"member operator{op}")
    else:  # free operator
        if op in _JL_OP and nargs == 2:
            parts, names = [], []
            for i, a in enumerate(entry["args"]):
                jt = jl_param_type(a["type"], wrapped)
                nm = a["name"] or f"a{i}"
                parts.append(f"{nm}::{jt}" if jt else nm)
                names.append(unwrap_arg(a["type"], nm, wrapped, nm in shifts))
            call = f"{op}({names[0]}, {names[1]})"
            expr, _ = wrap_return(ret, call, wrapped, ret_shift)
            return f"Base.:{_JL_OP[op]}({', '.join(parts)}) = {expr}"
        raise Skip(f"free operator{op}")


def emit_free_jl(entry, wrapped):
    name = entry["name"]
    ret_shift = entry["qualified"] in ov.ONEBASED_RETURN
    decls, callargs, names = build_params(entry, wrapped)
    call = f"{name}({', '.join(callargs)})"
    expr, _ = wrap_return(entry.get("return", {"cat": "void"}), call, wrapped,
                          ret_shift)
    # Pure passthrough (same name, no arg/return transform) would just recurse;
    # the CxxWrap function is directly usable, so no forwarder is needed.
    if callargs == names and expr == call:
        raise Skip("passthrough (use cxx directly)")
    return f"{name}(" + ", ".join(decls) + f") = {expr}"


def emit_field_jl(entry, wrapped):
    rec = short(entry["record"])
    fld = field_name(entry["qualified"].rsplit("::", 1)[0])
    name = entry["name"]
    call = f"{name}(r.{fld})"
    expr, _ = wrap_return(entry.get("type", {"cat": "void"}), call, wrapped,
                          False)
    return f"{name}(r::{rec}) = {expr}"


def emit_iterator(block, itr, wrapped):
    """Base.iterate for a block whose elements are ProductStates."""
    return f"""
function Base.iterate(b::{block})
    (size(b) > 0) || return nothing
    it = _begin(b.cxx_block)
    return ProductState(_deref(it)), it
end
function Base.iterate(b::{block}, it)
    _incr(it)
    it != _end(b.cxx_block) ? (ProductState(_deref(it)), it) : nothing
end"""


def _sig_key(line):
    """(funcname, (types...)) for a `NAME(params) = rhs` def.

    Two C++ overloads can collapse to the same Julia signature (e.g. the
    std::vector<int64_t> and arma::Col<int64_t> Permutation constructors both
    become Permutation(::Vector{Int64})). Julia would silently let the last
    definition win; here we dedup on the signature and keep the FIRST emitted
    one, which carries the correct index shift / return handling. Falls back to
    the whole line for anything that isn't a single `f(...) = ...` def."""
    head = line.split(" = ", 1)[0].strip()
    if not head.endswith(")"):
        return line
    depth, open_idx = 0, -1
    for i in range(len(head) - 1, -1, -1):
        if head[i] == ")":
            depth += 1
        elif head[i] == "(":
            depth -= 1
            if depth == 0:
                open_idx = i
                break
    if open_idx < 0:
        return line
    name = head[:open_idx]
    depth, cur, parts = 0, "", []
    for c in head[open_idx + 1:-1]:
        if c in "([{":
            depth += 1; cur += c
        elif c in ")]}":
            depth -= 1; cur += c
        elif c == "," and depth == 0:
            parts.append(cur); cur = ""
        else:
            cur += c
    if cur.strip():
        parts.append(cur)
    types = []
    for p in parts:
        p = p.split("=", 1)[0].strip()
        types.append(p.split("::", 1)[1].strip() if "::" in p else "Any")
    return (name, tuple(types))


def group_of(header):
    if header and header.startswith("xdiag/"):
        p = header.split("/")
        if len(p) >= 3:
            return p[1]
    return "misc"


def emit_struct(qual):
    name = short(qual)
    fld = field_name(qual)
    sup = " <: Block" if name in BLOCK_TYPES else ""
    return (f"struct {name}{sup}\n    {fld}::cxx_{name}\nend\n"
            f"convert(::Type{{T}}, x::cxx_{name}) where {{T<:{name}}} = {name}(x)")


def build(model):
    template_records = {e["qualified"] for e in model["entries"]
                        if e["kind"] == "record"
                        and e.get("record_kind") == "class_template"}
    wrapped = set()
    for q in model["api_record_names"]:
        if q in ov.EXCLUDE_RECORDS or q in ov.NATIVE_BRIDGES:
            continue
        if q in ov.RESULT_STRUCTS:
            continue  # materialised as native structs by the ergonomic layer
        if q in template_records or "Distributed" in q or short(q) in ITERATORS:
            continue
        wrapped.add(q)

    skipped = {}
    def note(q, why):
        skipped.setdefault(str(why), []).append(q)

    # All structs (+ abstract Block) go into one types file so forwarder files
    # can be included in any order.
    structs = [emit_struct(q) for q in sorted(wrapped, key=short)]

    methods = {}
    seen = set()
    def add(line, group):
        key = _sig_key(line)
        if key in seen:
            return
        seen.add(key)
        methods.setdefault(group, []).append(line)

    def wrappable_rec(rec):
        return rec and any(short(w) == rec for w in wrapped)

    for e in model["entries"]:
        if e["kind"] in ("record", "destructor"):
            continue
        q = e["qualified"]
        if q in ov.EXCLUDE_SYMBOLS or q in ov.SPECIALS or e.get("is_template"):
            note(q, "special/template/excluded")
            continue
        if e["kind"] in ("function", "method") and e.get("name") in ov.ERGONOMIC:
            note(q, "ergonomic (hand-written)")
            continue
        if q in ov.EXCLUDE_JULIA:
            note(q, "julia-excluded (hand-written kwargs forwarder)")
            continue
        g = group_of(e["header"])
        rec = e.get("record")
        # Free ops on native-bridge value types duplicate Julia's native ops
        # (and would recurse / clash with Base). Skip unless they yield a
        # wrapped type (e.g. scalar * Op -> OpSum).
        if not rec:
            _has_bridge = any(a["type"].get("cat") == "wrapped"
                              and a["type"].get("type") in ov.NATIVE_BRIDGES
                              for a in e.get("args", []))
            _ret = e.get("return", {"cat": "void"})
            _ret_wrapped = (_ret.get("cat") == "wrapped"
                            and _ret.get("type") not in ov.NATIVE_BRIDGES)
            if _has_bridge and not _ret_wrapped:
                note(q, "native-bridge value op")
                continue
        try:
            if e["kind"] == "constructor" and wrappable_rec(rec):
                add(emit_constructor(e, wrapped), g)
            elif e["kind"] == "method" and wrappable_rec(rec):
                add(emit_method(e, wrapped), g)
            elif e["kind"] == "var" and wrappable_rec(rec):
                add(emit_field_jl(e, wrapped), g)
            elif e["kind"] == "operator" and wrappable_rec(rec):
                add(emit_operator(e, wrapped, member=True), g)
            elif e["kind"] == "operator":
                add(emit_operator(e, wrapped, member=False), g)
            elif e["kind"] == "function":
                add(emit_free_jl(e, wrapped), g)
            else:
                note(q, e["kind"])
        except Skip as s:
            note(q, str(s))

    for e in model["entries"]:
        if e["kind"] == "record" and short(e["qualified"]) in BLOCK_TYPES:
            methods.setdefault(group_of(e["header"]), []).append(
                emit_iterator(short(e["qualified"]), None, wrapped))

    return structs, methods, skipped, wrapped


_HDR = ("# SPDX-License-Identifier: Apache-2.0\n"
        "# GENERATED by misc/wrapgen/emit_julia.py -- DO NOT EDIT.\n\n")


def render(model):
    structs, methods, skipped, wrapped = build(model)
    files = {}
    files["types.jl"] = (_HDR + "abstract type Block end\n\n"
                         + "\n\n".join(structs) + "\n")
    for g in sorted(methods):
        files[f"{g}.jl"] = _HDR + "\n".join(methods[g]) + "\n"
    return files, sorted(methods), skipped, wrapped


def main():
    import argparse, glob
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--report", action="store_true")
    args = ap.parse_args()
    model = _model.build_model([])
    errs = ov.validate(model)
    if errs:
        print("OVERRIDES INVALID:", *errs, sep="\n  ", file=sys.stderr); sys.exit(1)
    files, groups, skipped, wrapped = render(model)
    os.makedirs(args.out_dir, exist_ok=True)
    for old in glob.glob(os.path.join(args.out_dir, "*.jl")):
        os.remove(old)
    for fn, content in files.items():
        open(os.path.join(args.out_dir, fn), "w").write(content)
    print(f"Wrote {len(files)} Julia files ({len(wrapped)} types); groups: "
          + ", ".join(groups), file=sys.stderr)
    if args.report:
        tot = sum(len(v) for v in skipped.values())
        print(f"Skipped {tot}:", file=sys.stderr)
        for why, s in sorted(skipped.items(), key=lambda kv: -len(kv[1])):
            print(f"  {len(s):4d} {why}", file=sys.stderr)


if __name__ == "__main__":
    main()
