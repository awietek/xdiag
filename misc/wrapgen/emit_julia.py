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
            return name  # C++ overloads accept the native value directly
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
        return f"to_armadillo({name})"
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
            return expr, True  # already a native value from the C++ side
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


def build_params(entry, wrapped, self_type=None, self_field=None):
    """Return (param_decls, callargs). Raises Skip on an unmappable arg."""
    decls, callargs = [], []
    if self_type is not None:
        decls.append(f"self::{self_type}")
    shifts = onebased_args(entry["qualified"])
    for i, a in enumerate(entry.get("args", [])):
        jt = jl_param_type(a["type"], wrapped)
        nm = a["name"] or f"a{i}"
        decl = f"{nm}::{jt}" if jt else nm
        dv = default_jl(a.get("default")) if a.get("has_default") else None
        if dv is not None:
            decl += f"={dv}"
        decls.append(decl)
        callargs.append(unwrap_arg(a["type"], nm, wrapped, nm in shifts))
    return decls, callargs


def emit_constructor(entry, wrapped):
    rec = short(entry["record"])
    decls, callargs = build_params(entry, wrapped)
    call = f"construct_{rec}({', '.join(callargs)})"
    return f"{rec}({', '.join(decls)}) = {rec}({call})"


def emit_method(entry, wrapped):
    rec = short(entry["record"])
    fld = field_name(entry["qualified"].rsplit("::", 1)[0])
    ret_shift = entry["qualified"] in ov.ONEBASED_RETURN
    decls, callargs = build_params(entry, wrapped, self_type=rec)
    name = entry["name"]
    if entry["kind"] == "operator":
        raise Skip("operator (handled elsewhere)")
    call = f"{name}(self.{fld}{', ' if callargs else ''}{', '.join(callargs)})"
    wrapped_expr, _ = wrap_return(entry.get("return", {"cat": "void"}), call,
                                  wrapped, ret_shift)
    return f"{name}(" + ", ".join(decls) + f") = {wrapped_expr}"


def group_of(header):
    if header and header.startswith("xdiag/"):
        parts = header.split("/")
        if len(parts) >= 3:
            return parts[1]
    return "misc"


def emit_struct(qual):
    """Julia wrapper struct + convert for a wrapped type."""
    name = short(qual)
    fld = field_name(qual)
    return (f"struct {name}\n    {fld}::cxx_{name}\nend\n"
            f"convert(::Type{{T}}, x::cxx_{name}) where {{T<:{name}}} = {name}(x)")


def emit_iterator(block, itr, wrapped):
    """Base.iterate protocol for a block whose elements are ProductStates."""
    fld = field_name(next(q for q in wrapped if short(q) == block))
    return f"""
function Base.iterate(b::{block})
    if size(b) > 0
        it = begin(b.{fld})
        return ProductState(_deref(it)), it
    end
    return nothing
end
function Base.iterate(b::{block}, it)
    _incr(it)
    it != end(b.{fld}) ? (ProductState(_deref(it)), it) : nothing
end"""


def build(model):
    template_records = {e["qualified"] for e in model["entries"]
                        if e["kind"] == "record"
                        and e.get("record_kind") == "class_template"}
    wrapped = set()
    for q in model["api_record_names"]:
        if q in ov.EXCLUDE_RECORDS or q in ov.NATIVE_BRIDGES:
            continue
        if q in template_records or "Distributed" in q or short(q) in ITERATORS:
            continue
        wrapped.add(q)

    skipped = {}
    def note(q, why):
        skipped.setdefault(str(why), []).append(q)

    # structs (all wrapped types, grouped by subsystem)
    struct_group, method_group = {}, {}
    rec_header = {}
    for e in model["entries"]:
        if e["kind"] == "record" and e["qualified"] in wrapped:
            g = group_of(e["header"])
            rec_header[short(e["qualified"])] = g
            struct_group.setdefault(g, []).append(emit_struct(e["qualified"]))

    seen = set()
    def add(line, name, group):
        key = (group, line.split("(")[0], line)
        if key in seen:
            return
        seen.add(key)
        method_group.setdefault(group, []).append(line)

    for e in model["entries"]:
        if e["kind"] == "record" or e["kind"] == "var":
            continue
        q = e["qualified"]
        if q in ov.EXCLUDE_SYMBOLS or q in ov.SPECIALS or e.get("is_template"):
            note(q, "special/template/excluded")
            continue
        rec = e.get("record")
        g = group_of(e["header"])
        try:
            if e["kind"] == "constructor":
                if rec and any(short(w) == rec for w in wrapped):
                    add(emit_constructor(e, wrapped), rec, g)
            elif e["kind"] == "method":
                if rec and any(short(w) == rec for w in wrapped):
                    add(emit_method(e, wrapped), rec, g)
            else:
                note(q, e["kind"])
        except Skip as s:
            note(q, str(s))

    # iterators
    for e in model["entries"]:
        if e["kind"] == "record" and short(e["qualified"]) in BLOCK_TYPES:
            g = group_of(e["header"])
            method_group.setdefault(g, []).append(
                emit_iterator(short(e["qualified"]), None, wrapped))

    return struct_group, method_group, skipped, wrapped


_FILE_HDR = "# SPDX-License-Identifier: Apache-2.0\n# GENERATED by misc/wrapgen/emit_julia.py -- DO NOT EDIT.\n\n"


def render(model):
    structs, methods, skipped, wrapped = build(model)
    files = {}
    groups = sorted(set(structs) | set(methods))
    for g in groups:
        body = _FILE_HDR
        body += "\n\n".join(structs.get(g, [])) + "\n\n"
        body += "\n".join(methods.get(g, [])) + "\n"
        files[f"{g}.jl"] = body
    return files, groups, skipped, wrapped


def main():
    import argparse
    import glob
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--report", action="store_true")
    args = ap.parse_args()

    model = _model.build_model([])
    errs = ov.validate(model)
    if errs:
        print("OVERRIDES INVALID:", *errs, sep="\n  ", file=sys.stderr)
        sys.exit(1)

    files, groups, skipped, wrapped = render(model)
    os.makedirs(args.out_dir, exist_ok=True)
    for old in glob.glob(os.path.join(args.out_dir, "*.jl")):
        os.remove(old)
    for fn, content in files.items():
        with open(os.path.join(args.out_dir, fn), "w") as fh:
            fh.write(content)
    print(f"Wrote {len(files)} Julia files to {args.out_dir} "
          f"({len(wrapped)} types).", file=sys.stderr)
    if args.report:
        total = sum(len(v) for v in skipped.values())
        print(f"Skipped {total}:", file=sys.stderr)
        for why, s in sorted(skipped.items(), key=lambda kv: -len(kv[1])):
            print(f"  {len(s):4d} {why}", file=sys.stderr)


if __name__ == "__main__":
    main()
