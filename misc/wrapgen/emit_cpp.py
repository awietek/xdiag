# SPDX-License-Identifier: Apache-2.0
#
# Emit the CxxWrap C++ side of the Julia wrapper from the API model + overrides.
#
# Design: the C++ side is a thin pass-through. Every non-excluded XDIAG_API
# symbol becomes a `.method(...)` / `mod.method(...)` registration whose lambda
# just calls the underlying C++ entity (wrapped in the JULIA_XDIAG_CALL_* macros
# for stack-trace-preserving error translation). All 1-based<->0-based index
# shifts, native-type bridging and default arguments are handled on the Julia
# side (emit_julia.py), so they do not appear here.
#
# Registration is TWO-PASS: every add_type<> is emitted before any method, which
# removes the "No appropriate factory" ordering fragility entirely.
#
# This is a first cut: it emits the mechanical majority (records + their
# methods/constructors/common operators, and free functions) whose args/returns
# are all trivially marshallable. Symbols needing hand-written specials, native
# bridges, template instantiation, or categories not yet handled are skipped and
# reported so coverage stays honest.

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import model as _model
import overrides as ov

# Canonical arma spelling -> the arma typedef used in lambdas.
_ARMA_CPP = {
    "vec": "arma::vec", "cx_vec": "arma::cx_vec",
    "mat": "arma::mat", "cx_mat": "arma::cx_mat",
    "ivec": "arma::Col<int64_t>", "imat": "arma::Mat<int64_t>",
}

# Primitive canonical spelling -> C++ lambda type.
_PRIM_CPP = {
    "bool": "bool", "int32": "int32_t", "int64": "int64_t",
    "uint32": "uint32_t", "uint64": "uint64_t", "int16": "int16_t",
    "double": "double", "float": "float", "complex": "complex", "char": "char",
}

# Binary operators we emit (Julia method name == the symbol). Mutating and
# ordering operators are skipped in this first cut.
_BINARY_OPS = {"+", "-", "*", "/", "==", "!="}

# Concrete serial block types the polymorphic Block variant is expanded over.
BLOCK_TYPES = ["Spinhalf", "tJ", "Electron", "Boson", "Fermion"]

# Plain aggregate result structs CxxWrap treats as mirrored types: they cannot
# be add_type'd and should instead be mapped to native Julia structs (later).
MIRRORED = {
    "xdiag::EigsLanczosResult", "xdiag::EigvalsLanczosResult",
    "xdiag::EvolveLanczosResult", "xdiag::EvolveLanczosInplaceResult",
    "xdiag::TimeEvolveExpokitResult", "xdiag::TimeEvolveExpokitInplaceResult",
    "xdiag::EigsLobpcgResult",
}


def short(qualified):
    return qualified.split("::")[-1]


class Skip(Exception):
    pass


def cpp_type(tp, wrapped_names, block_sub=None):
    """C++ lambda parameter/return type for a classified type, or raise Skip if
    this category is not handled by the mechanical emitter yet. block_sub is the
    concrete block name substituted for a Block-variant parameter."""
    cat = tp["cat"]
    ref = tp.get("ref")
    const = tp.get("const", False)

    def decorate(base, value_ok=True):
        # Reproduce const/ref of the original parameter where it matters.
        if ref == "lvalue":
            return f"{base} const &" if const else f"{base} &"
        if ref == "rvalue":
            return f"{base} &&"
        return base

    if cat == "void":
        return "void"
    if cat == "primitive":
        return decorate(_PRIM_CPP.get(tp.get("prim"), None) or raise_skip(tp))
    if cat == "string":
        return decorate("std::string")
    if cat == "arma":
        return decorate(_ARMA_CPP.get(tp.get("arma")) or raise_skip(tp))
    if cat == "wrapped":
        name = tp["type"]
        if name in ov.EXCLUDE_RECORDS or name in ov.NATIVE_BRIDGES:
            raise Skip(f"wrapped/bridge type {name}")
        if name not in wrapped_names:
            raise Skip(f"unwrapped type {name}")
        return decorate(short(name))
    if cat == "variant":  # the polymorphic Block
        if tp.get("is_block"):
            if block_sub is None:
                raise Skip("block variant (no substitution)")
            return decorate(block_sub)
        raise Skip("non-block variant")
    if cat == "stdvector":
        elem = tp.get("elem", {})
        ecat = elem.get("cat")
        if ecat == "primitive":
            return decorate(f"std::vector<{_PRIM_CPP.get(elem.get('prim'),'int64_t')}>")
        if ecat == "string":
            return decorate("std::vector<std::string>")
        if ecat == "wrapped" and elem["type"] in wrapped_names:
            return decorate(f"std::vector<{short(elem['type'])}>")
        raise Skip(f"stdvector<{ecat}>")
    if cat == "enum":
        raise Skip("enum")
    raise Skip(cat)


def raise_skip(tp):
    raise Skip(tp.get("cat", "?"))


def has_block_param(entry):
    return any(a["type"].get("cat") == "variant" and a["type"].get("is_block")
               for a in entry.get("args", []))


def block_subs(entry):
    """Concrete block names to expand over ([None] if no Block params)."""
    return BLOCK_TYPES if has_block_param(entry) else [None]


def lambda_params(entry, wrapped_names, self_type=None, block_sub=None):
    """Return (param_decls, arg_names). Prepends a self param for members."""
    decls, names = [], []
    if self_type is not None:
        const = entry.get("is_const", True)
        decls.append(f"{self_type} {'const &' if const else '&'} self")
        names.append("self")
    for i, a in enumerate(entry.get("args", [])):
        t = cpp_type(a["type"], wrapped_names, block_sub)
        if t == "void":
            raise Skip("void param")
        an = a["name"] or f"a{i}"
        decls.append(f"{t} {an}")
        names.append(an)
    return decls, names


def call_wrap(entry, expr):
    ret = entry.get("return", {"cat": "void"})
    if ret["cat"] == "void":
        return f"JULIA_XDIAG_CALL_VOID({expr});"
    return f"JULIA_XDIAG_CALL_RETURN({expr})"


def _check_return(entry, wrapped_names):
    """Validate the return type marshals; variant(block) returns are skipped."""
    ret = entry.get("return", {"cat": "void"})
    if ret.get("cat") == "variant":
        raise Skip("returns Block variant")
    cpp_type(ret, wrapped_names)  # raises Skip if not marshallable


def emit_method(tw, entry, wrapped_names):
    """A member method / operator on a wrapped record. Returns a list of lines
    (one per Block substitution)."""
    rec = short(entry["qualified"].rsplit("::", 1)[0])
    args = entry.get("args", [])
    _check_return(entry, wrapped_names)
    out = []
    for sub in block_subs(entry):
        decls, names = lambda_params(entry, wrapped_names, self_type=rec,
                                     block_sub=sub)
        call_args = ", ".join(names[1:])
        if entry["kind"] == "operator":
            op = entry.get("op")
            if op in _BINARY_OPS and len(args) == 1:
                jl_name, expr = op, f"self {op} {names[1]}"
            elif op == "[]" and len(args) == 1:
                jl_name, expr = "getindex", f"self[{names[1]}]"
            else:
                raise Skip(f"operator{op}")
        else:
            jl_name = entry["name"]
            expr = f"self.{entry['name']}({call_args})"
        body = call_wrap(entry, expr)
        out.append(f'  {tw}.method("{jl_name}", []({", ".join(decls)}) '
                   f'{{ {body} }});')
    return out


def emit_constructor(entry, wrapped_names):
    rec = short(entry["qualified"].rsplit("::", 1)[0])
    out = []
    for sub in block_subs(entry):
        decls, names = lambda_params(entry, wrapped_names, block_sub=sub)
        out.append(f'  mod.method("construct_{rec}", []({", ".join(decls)}) '
                   f'{{ JULIA_XDIAG_CALL_RETURN({rec}({", ".join(names)})) }});')
    return out


def emit_free(entry, wrapped_names):
    name = entry["name"]
    _check_return(entry, wrapped_names)
    out = []
    for sub in block_subs(entry):
        decls, names = lambda_params(entry, wrapped_names, block_sub=sub)
        if entry["kind"] == "operator":
            op = entry.get("op")
            if op in _BINARY_OPS and len(names) == 2:
                jl_name, expr = op, f"{names[0]} {op} {names[1]}"
            else:
                raise Skip(f"free operator{op}")
        else:
            jl_name = name
            expr = f"{name}({', '.join(names)})"
        body = call_wrap(entry, expr)
        out.append(f'  mod.method("{jl_name}", []({", ".join(decls)}) '
                   f'{{ {body} }});')
    return out


def excluded(entry):
    header = entry.get("header") or ""
    for pat in ov.EXCLUDE_HEADERS:
        if pat.replace("*", "") in header:
            return True
    q = entry.get("qualified", "")
    if q in ov.EXCLUDE_SYMBOLS:
        return True
    if q in ov.SPECIALS:
        return True  # emitted from a hand-written fragment instead
    return False


def build(model):
    # Class-template records (CSRMatrix<>, ...) cannot be add_type'd generically;
    # specific instantiations are provided by hand-written specials.
    template_records = {e["qualified"] for e in model["entries"]
                        if e["kind"] == "record"
                        and e.get("record_kind") == "class_template"}

    wrapped_names = set()
    for q in model["api_record_names"]:
        if q in ov.EXCLUDE_RECORDS or q in ov.NATIVE_BRIDGES or q in MIRRORED:
            continue
        if q in template_records or "Distributed" in q:
            continue
        wrapped_names.add(q)

    tw_of = {}  # qualified record name -> TypeWrapper variable
    lines_pass1, lines_pass2 = [], []
    skipped = {}

    def note_skip(q, why):
        skipped.setdefault(str(why), []).append(q)

    # PASS 1: register every wrapped type.
    for q in sorted(wrapped_names):
        name = short(q)
        field = ov.FIELD_NAMES.get(q, "cxx_" + _snake(name))
        var = "tw_" + _snake(name)
        tw_of[q] = var
        lines_pass1.append(f'  auto {var} = mod.add_type<{name}>("cxx_{name}");')

    # PASS 2: methods, constructors, operators, free functions.
    records = {}
    frees = []
    for e in model["entries"]:
        if e["kind"] == "record":
            continue
        if excluded(e):
            continue
        rec = e.get("record")
        if rec:
            recq = next((q for q in wrapped_names if short(q) == rec), None)
            if recq is None:
                note_skip(e["qualified"], "record not wrapped")
                continue
            records.setdefault(recq, []).append(e)
        elif e["kind"] in ("function", "operator"):
            frees.append(e)

    for q in sorted(records):
        tw = tw_of[q]
        lines_pass2.append(f"  // ---- {short(q)} ----")
        for e in sorted(records[q], key=lambda e: e["line"]):
            try:
                if e["kind"] == "constructor":
                    lines_pass2.extend(emit_constructor(e, wrapped_names))
                elif e["kind"] in ("method", "operator"):
                    lines_pass2.extend(emit_method(tw, e, wrapped_names))
                elif e["kind"] == "destructor":
                    pass
                else:
                    note_skip(e["qualified"], "member kind " + e["kind"])
            except Skip as s:
                note_skip(e["qualified"], str(s))

    lines_pass2.append("  // ---- free functions ----")
    for e in sorted(frees, key=lambda e: (e["header"] or "", e["line"])):
        try:
            lines_pass2.extend(emit_free(e, wrapped_names))
        except Skip as s:
            note_skip(e["qualified"], str(s))

    return lines_pass1, lines_pass2, skipped, wrapped_names


def _snake(name):
    out = []
    for i, c in enumerate(name):
        if c.isupper() and i > 0 and not name[i - 1].isupper():
            out.append("_")
        out.append(c.lower())
    return "".join(out)


HEADER = """// SPDX-License-Identifier: Apache-2.0
// GENERATED by misc/wrapgen/emit_cpp.py -- DO NOT EDIT.
// Source model: XDIAG_API surface. Regenerate with misc/wrapgen/generate.py.

#include <julia/src/xdiagjl.hpp>

namespace xdiag::julia {

void define_generated(jlcxx::Module &mod) {
  using namespace xdiag;

"""

FOOTER = """}

} // namespace xdiag::julia
"""


def render(model):
    p1, p2, skipped, wrapped = build(model)
    body = HEADER
    body += "  // ==== PASS 1: register all types ====\n"
    body += "\n".join(p1) + "\n\n"
    body += "  // ==== PASS 2: methods / constructors / free functions ====\n"
    body += "\n".join(p2) + "\n"
    body += FOOTER
    return body, skipped, wrapped


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default="-")
    ap.add_argument("--report", action="store_true")
    args = ap.parse_args()

    model = _model.build_model([])
    errs = ov.validate(model)
    if errs:
        print("OVERRIDES INVALID:", *errs, sep="\n  ", file=sys.stderr)
        sys.exit(1)

    code, skipped, wrapped = render(model)
    if args.output == "-":
        sys.stdout.write(code)
    else:
        with open(args.output, "w") as f:
            f.write(code)
        print(f"Wrote {args.output} ({len(wrapped)} types, "
              f"{code.count('.method(') + code.count('mod.method(')} registrations).",
              file=sys.stderr)

    if args.report:
        total = sum(len(v) for v in skipped.values())
        print(f"\nSkipped {total} symbols by reason:", file=sys.stderr)
        for why, syms in sorted(skipped.items(), key=lambda kv: -len(kv[1])):
            print(f"  {len(syms):4d}  {why}", file=sys.stderr)


if __name__ == "__main__":
    main()
