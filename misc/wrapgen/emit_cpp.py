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


def emit_method(entry, wrapped_names):
    """Member method / operator on a wrapped record -> list of (sig_key, line).
    Registered via mod.method (self is the explicit first arg) so it needs no
    TypeWrapper -- this lets PASS 2 be split across files."""
    rec = short(entry["qualified"].rsplit("::", 1)[0])
    nargs = len(entry.get("args", []))
    wrap, is_void = return_transform(entry, wrapped_names)
    out = []
    for decls, callargs, sig in expansions(entry, wrapped_names, self_type=rec):
        void = is_void
        if entry["kind"] == "operator":
            op = entry.get("op")
            if op in _BINARY_OPS and nargs == 1:
                jl_name, expr = op, f"self {op} {callargs[0]}"
            elif op == "[]" and nargs == 1:
                jl_name, expr = "getindex", f"self[{callargs[0]}]"
            elif op == "++" and nargs == 0:  # iterator increment
                jl_name, expr, void = "_incr", "++self", True
            elif op == "*" and nargs == 0:  # iterator dereference
                jl_name, expr = "_deref", "*self"
            else:
                raise Skip(f"operator{op}")
        else:
            jl_name = entry["name"]
            expr = f"self.{entry['name']}({', '.join(callargs)})"
        expr = wrap(expr)
        body = (f"JULIA_XDIAG_CALL_VOID({expr});" if void
                else f"JULIA_XDIAG_CALL_RETURN({expr})")
        out.append(((jl_name, sig), _reg("mod", jl_name, decls, body)))
    return out


def emit_constructor(entry, wrapped_names):
    rec = short(entry["qualified"].rsplit("::", 1)[0])
    out = []
    for decls, callargs, sig in expansions(entry, wrapped_names):
        body = f'JULIA_XDIAG_CALL_RETURN({rec}({", ".join(callargs)}))'
        out.append(((f"construct_{rec}", sig),
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
        out.append(((jl_name, sig), _reg("mod", jl_name, decls, body)))
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
    # Per-overload exclusion: drop this overload if an argument type matches.
    subs = ov.EXCLUDE_OVERLOADS.get(q)
    if subs:
        for a in entry.get("args", []):
            t = a["type"]
            s = t.get("type") or t.get("cpp") or ""
            if any(sub in s for sub in subs):
                return True
    return False


def _snake(name):
    out = []
    for i, c in enumerate(name):
        if c.isupper() and i > 0 and not name[i - 1].isupper():
            out.append("_")
        out.append(c.lower())
    return "".join(out)


def group_of(header):
    """Subsystem bucket for a symbol, from its header path (xdiag/<group>/...)."""
    if header and header.startswith("xdiag/"):
        parts = header.split("/")
        if len(parts) >= 3:
            return parts[1]
    return "misc"


def build(model):
    # Class-template records (CSRMatrix<>, ...) can't be add_type'd generically.
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

    skipped = {}

    def note_skip(q, why):
        skipped.setdefault(str(why), []).append(q)

    # PASS 1: register every wrapped type (no TypeWrapper kept -- all methods go
    # through mod.method, so PASS 2 can live in separate translation units).
    types = [f'  mod.add_type<{short(q)}>("cxx_{short(q)}");'
             for q in sorted(wrapped_names)]

    # PASS 2: bucket registrations by subsystem so they compile in parallel.
    records, frees = {}, []
    for e in model["entries"]:
        if e["kind"] == "record" or excluded(e):
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

    seen = set()
    groups = {}

    def add(regs, q, group):
        for key, line in regs:
            jl_name, sig = key
            nkey = (jl_name, tuple("".join(s.split()) for s in sig))
            if nkey in seen:
                note_skip(q, "duplicate signature")
                continue
            seen.add(nkey)
            groups.setdefault(group, []).append(line)

    for q in sorted(records):
        for e in sorted(records[q], key=lambda e: e["line"]):
            group = group_of(e["header"])
            try:
                if e["kind"] == "constructor":
                    add(emit_constructor(e, wrapped_names), e["qualified"], group)
                elif e["kind"] in ("method", "operator"):
                    add(emit_method(e, wrapped_names), e["qualified"], group)
                elif e["kind"] == "destructor":
                    pass
                else:
                    note_skip(e["qualified"], "member kind " + e["kind"])
            except Skip as s:
                note_skip(e["qualified"], str(s))

    for e in sorted(frees, key=lambda e: (e["header"] or "", e["line"])):
        try:
            add(emit_free(e, wrapped_names), e["qualified"], group_of(e["header"]))
        except Skip as s:
            note_skip(e["qualified"], str(s))

    return types, groups, skipped, wrapped_names


_FILE_HEADER = """// SPDX-License-Identifier: Apache-2.0
// GENERATED by misc/wrapgen/emit_cpp.py -- DO NOT EDIT.

#include <julia/src/xdiagjl.hpp>

namespace xdiag::julia {
"""
_FILE_FOOTER = "\n} // namespace xdiag::julia\n"


def _fn(name, body_lines):
    return (f"\nvoid {name}(jlcxx::Module &mod) {{\n"
            f"  using namespace xdiag;\n" + "\n".join(body_lines) + "\n}\n")


def render(model):
    types, groups, skipped, wrapped = build(model)
    files = {}
    files["define_types.cpp"] = (_FILE_HEADER + _fn("define_types", types)
                                 + _FILE_FOOTER)
    gnames = sorted(groups)
    for g in gnames:
        files[f"define_{g}.cpp"] = (_FILE_HEADER + _fn(f"define_{g}", groups[g])
                                    + _FILE_FOOTER)
    decls = ["void define_types(jlcxx::Module &mod);"]
    calls = ["  define_types(mod);"]
    for g in gnames:
        decls.append(f"void define_{g}(jlcxx::Module &mod);")
        calls.append(f"  define_{g}(mod);")
    files["module_generated.cpp"] = (
        _FILE_HEADER + "\n" + "\n".join(decls)
        + "\n\nvoid define_generated(jlcxx::Module &mod) {\n"
        + "\n".join(calls) + "\n}\n" + _FILE_FOOTER)
    return files, skipped, wrapped


def main():
    import argparse
    import glob
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", help="directory for the generated .cpp files")
    ap.add_argument("--report", action="store_true")
    args = ap.parse_args()

    model = _model.build_model([])
    errs = ov.validate(model)
    if errs:
        print("OVERRIDES INVALID:", *errs, sep="\n  ", file=sys.stderr)
        sys.exit(1)

    files, skipped, wrapped = render(model)

    if args.out_dir:
        os.makedirs(args.out_dir, exist_ok=True)
        # Clear stale generated files first (a subsystem may have vanished).
        for old in glob.glob(os.path.join(args.out_dir, "define_*.cpp")):
            os.remove(old)
        stale = os.path.join(args.out_dir, "module_generated.cpp")
        if os.path.exists(stale):
            os.remove(stale)
        nreg = 0
        for fn, content in files.items():
            with open(os.path.join(args.out_dir, fn), "w") as fh:
                fh.write(content)
            nreg += content.count("mod.method(")
        print(f"Wrote {len(files)} files to {args.out_dir} "
              f"({len(wrapped)} types, {nreg} registrations).", file=sys.stderr)
    else:
        for fn, content in files.items():
            sys.stdout.write(f"// ==== {fn} ====\n{content}\n")

    if args.report:
        total = sum(len(v) for v in skipped.values())
        print(f"\nSkipped {total} symbols by reason:", file=sys.stderr)
        for why, syms in sorted(skipped.items(), key=lambda kv: -len(kv[1])):
            print(f"  {len(syms):4d}  {why}", file=sys.stderr)


if __name__ == "__main__":
    main()
