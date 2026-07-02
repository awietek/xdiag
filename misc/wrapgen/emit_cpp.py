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

import itertools
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

# Native-type bridges: a bridge parameter is accepted as each of these C++
# types (by value) and wrapped into the bridge type at the call site; a bridge
# return is converted back. The Julia side then sees plain Number/String/arrays.
NATIVE_BRIDGE_PARAM = {
    "scalar": [("double", "Scalar({n})"), ("complex", "Scalar({n})")],
    "coeff": [("std::string", "Coeff({n})"), ("double", "Coeff({n})"),
              ("complex", "Coeff({n})")],
    "vector": [("arma::vec", "Vector({n})"), ("arma::cx_vec", "Vector({n})")],
    "matrix": [("arma::mat", "Matrix({n})"), ("arma::cx_mat", "Matrix({n})")],
}
NATIVE_BRIDGE_RETURN = {  # expr -> converted expr; missing kinds are skipped
    "scalar": "({expr}).as<complex>()",
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


def param_choices(tp, wrapped_names):
    """Concrete lambda-parameter choices for one API parameter, as a list of
    (cpp_decl_type, argfmt) where argfmt.format(n=name) is the value passed to
    the C++ callee. Multiple entries mean the parameter is expanded: a Block
    over the concrete serial blocks, a native-bridge type over its accepted
    Julia representations."""
    cat = tp["cat"]
    if cat == "variant" and tp.get("is_block"):
        return [(f"{bt} const &", "{n}") for bt in BLOCK_TYPES]
    if cat == "wrapped" and tp["type"] in ov.NATIVE_BRIDGES:
        opts = NATIVE_BRIDGE_PARAM.get(ov.NATIVE_BRIDGES[tp["type"]])
        if not opts:
            raise Skip(f"bridge param {tp['type']}")
        return list(opts)
    return [(cpp_type(tp, wrapped_names), "{n}")]


def return_transform(entry, wrapped_names):
    """(expr_wrapper, is_void). Validates the return marshals or raises Skip.
    A native-bridge return is converted back to its Julia representation."""
    ret = entry.get("return", {"cat": "void"})
    cat = ret.get("cat")
    if cat == "void":
        return (lambda e: e), True
    if cat == "variant":
        raise Skip("returns Block variant")
    if cat == "wrapped" and ret["type"] in ov.NATIVE_BRIDGES:
        tmpl = NATIVE_BRIDGE_RETURN.get(ov.NATIVE_BRIDGES[ret["type"]])
        if not tmpl:
            raise Skip(f"bridge return {ret['type']}")
        return (lambda e: tmpl.format(expr=e)), False
    cpp_type(ret, wrapped_names)  # validate marshallable
    return (lambda e: e), False


def expansions(entry, wrapped_names, self_type=None):
    """Yield (decls, callargs, sig) for each concrete overload of a callable.
    sig is the tuple of C++ parameter types used for de-duplication."""
    choices = [param_choices(a["type"], wrapped_names)
               for a in entry.get("args", [])]
    for combo in itertools.product(*choices):
        decls, callargs, sig = [], [], []
        if self_type is not None:
            const = entry.get("is_const", True)
            decls.append(f"{self_type} {'const &' if const else '&'} self")
            sig.append(self_type + (" const&" if const else "&"))
        for i, (a, (decl, fmt)) in enumerate(zip(entry.get("args", []), combo)):
            an = a["name"] or f"a{i}"
            decls.append(f"{decl} {an}")
            callargs.append(fmt.format(n=an))
            sig.append(decl)
        yield decls, callargs, tuple(sig)


def _reg(target, jl_name, decls, body):
    return f'  {target}.method("{jl_name}", []({", ".join(decls)}) {{ {body} }});'


def emit_method(tw, entry, wrapped_names):
    """Member method / operator on a wrapped record -> list of (sig_key, line)."""
    rec = short(entry["qualified"].rsplit("::", 1)[0])
    nargs = len(entry.get("args", []))
    wrap, is_void = return_transform(entry, wrapped_names)
    out = []
    for decls, callargs, sig in expansions(entry, wrapped_names, self_type=rec):
        if entry["kind"] == "operator":
            op = entry.get("op")
            if op in _BINARY_OPS and nargs == 1:
                jl_name, expr = op, f"self {op} {callargs[0]}"
            elif op == "[]" and nargs == 1:
                jl_name, expr = "getindex", f"self[{callargs[0]}]"
            else:
                raise Skip(f"operator{op}")
        else:
            jl_name = entry["name"]
            expr = f"self.{entry['name']}({', '.join(callargs)})"
        expr = wrap(expr)
        body = (f"JULIA_XDIAG_CALL_VOID({expr});" if is_void
                else f"JULIA_XDIAG_CALL_RETURN({expr})")
        out.append(((tw, jl_name, sig), _reg(tw, jl_name, decls, body)))
    return out


def emit_constructor(entry, wrapped_names):
    rec = short(entry["qualified"].rsplit("::", 1)[0])
    out = []
    for decls, callargs, sig in expansions(entry, wrapped_names):
        body = f'JULIA_XDIAG_CALL_RETURN({rec}({", ".join(callargs)}))'
        out.append((("mod", f"construct_{rec}", sig),
                    _reg("mod", f"construct_{rec}", decls, body)))
    return out


def emit_free(entry, wrapped_names):
    name = entry["name"]
    nargs = len(entry.get("args", []))
    wrap, is_void = return_transform(entry, wrapped_names)
    out = []
    for decls, callargs, sig in expansions(entry, wrapped_names):
        if entry["kind"] == "operator":
            op = entry.get("op")
            if op in _BINARY_OPS and nargs == 2:
                jl_name, expr = op, f"{callargs[0]} {op} {callargs[1]}"
            else:
                raise Skip(f"free operator{op}")
        else:
            jl_name = name
            expr = f"{name}({', '.join(callargs)})"
        expr = wrap(expr)
        body = (f"JULIA_XDIAG_CALL_VOID({expr});" if is_void
                else f"JULIA_XDIAG_CALL_RETURN({expr})")
        out.append((("mod", jl_name, sig), _reg("mod", jl_name, decls, body)))
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

    seen = set()  # de-dup registrations whose C++ signatures collide

    def add(regs, q):
        for sig, line in regs:
            if sig in seen:
                note_skip(q, "duplicate signature")
                continue
            seen.add(sig)
            lines_pass2.append(line)

    for q in sorted(records):
        tw = tw_of[q]
        lines_pass2.append(f"  // ---- {short(q)} ----")
        for e in sorted(records[q], key=lambda e: e["line"]):
            try:
                if e["kind"] == "constructor":
                    add(emit_constructor(e, wrapped_names), e["qualified"])
                elif e["kind"] in ("method", "operator"):
                    add(emit_method(tw, e, wrapped_names), e["qualified"])
                elif e["kind"] == "destructor":
                    pass
                else:
                    note_skip(e["qualified"], "member kind " + e["kind"])
            except Skip as s:
                note_skip(e["qualified"], str(s))

    lines_pass2.append("  // ---- free functions ----")
    for e in sorted(frees, key=lambda e: (e["header"] or "", e["line"])):
        try:
            add(emit_free(e, wrapped_names), e["qualified"])
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
