#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
#
# Emit the callables per module (out/cpp/<module>.cpp, out/julia/<module>.jl):
# constructors, methods, free functions, the common operators, and struct
# field accessors. One binding per symbol -- the ABI name (con_/mth_/fun_/op_)
# is chosen once and used on both sides. Non-conforming symbols are skipped and
# listed in out/coverage.txt.

from collections import defaultdict

import apidb
import emitutil
import overrides
import typetable as T
from emitutil import register_method

BINARY_OPS = {"operator+": "+", "operator-": "-", "operator*": "*",
              "operator/": "/", "operator==": "=="}
OP_ABI = {"+": "add", "-": "sub", "*": "mul", "/": "div", "==": "eq"}
OP_JULIA = {"+": "Base.:+", "-": "Base.:-", "*": "Base.:*", "/": "Base.:/",
            "==": "Base.:(==)"}
BASE = {"size", "length", "isreal", "show", "real", "imag", "conj", "abs",
        "adjoint", "isapprox", "inv", "zero", "hash", "isvalid"}
LINALG = {"norm", "dot"}


def julia_name(name):
    if name in BASE:
        return f"Base.{name}"
    if name in LINALG:
        return f"LinearAlgebra.{name}"
    return name


def julia_default(cpp):
    if cpp is None:
        return None
    return {"true": "true", "false": "false", "nullptr": "nothing"}.get(cpp, cpp.rstrip("Ll"))


def cpp_self(e):
    return f"{e['record']} const& self" if e.get("const") else f"{e['record']}& self"


# concrete (idx, coeff) instantiations for a CSRMatrix<...> argument
SPARSE_COMBOS = [("int64_t", "double"), ("int64_t", "complex"),
                 ("int32_t", "double"), ("int32_t", "complex")]


def cpp_variants(args):
    """C++ parameter lists. A Block argument expands over the concrete blocks
    and a CSRMatrix argument over the (idx, coeff) instantiations."""
    blocks = T.BLOCKS if any(T.category(a["type"]) == "block" for a in args) else [None]
    sparses = SPARSE_COMBOS if any(T.category(a["type"]) == "sparse" for a in args) else [None]
    for block in blocks:
        for sparse in sparses:
            yield [T.cpp_param(a["type"], a["name"], block_sub=block, sparse_sub=sparse)
                   for a in args]


def julia_params(args, shifts=(), lead=None):
    """Julia signature string. C++ default arguments (always a trailing suffix)
    become keyword arguments after ';'. A defaulted index arg (in `shifts`) has
    its default bumped by 1 (C++ 0-based -> Julia 1-based)."""
    required = [lead] if lead else []
    optional = []
    for a in args:
        decl = f"{a['name']}::{T.julia_type(a['type'])}"
        if a["default"]:
            d = julia_default(a["default"])
            if a["name"] in shifts and d.lstrip("-").isdigit():
                d = str(int(d) + 1)
            optional.append(f"{decl} = {d}")
        else:
            required.append(decl)
    sig = ", ".join(required)
    if optional:
        sig += "; " + ", ".join(optional)
    return sig


def julia_call_args(args, shifts, lead=None):
    parts = [lead] if lead else []
    for a in args:
        value = T.shift_down(a["name"]) if a["name"] in shifts else a["name"]
        parts.append(T.unwrap(a["type"], value))
    return ", ".join(parts)


# --- one binding per symbol -------------------------------------------------
def bind_constructor(e, cpp, jl):
    cls, names = e["record"], ", ".join(a["name"] for a in e["args"])
    for params in cpp_variants(e["args"]):
        cpp.append(register_method(f"con_{cls}", params, f"{cls}({names})", False))
    jl.append(f"{cls}({julia_params(e['args'], overrides.arg_shifts(e))}) = "
              f"{cls}(con_{cls}({julia_call_args(e['args'], overrides.arg_shifts(e))}))")


def bind_method(e, cpp, jl):
    name, ret = e["name"], e.get("returns")
    if name in MERGED_NAMES:        # a merged real/complex free function supersedes it
        return
    names = ", ".join(a["name"] for a in e["args"])
    for params in cpp_variants(e["args"]):
        cpp.append(register_method(f"mth_{name}", [cpp_self(e)] + params,
                                   f"self.{name}({names})", ret is None))
    call = f"mth_{name}({julia_call_args(e['args'], overrides.arg_shifts(e), 'obj.cxx_object')})"
    result = T.wrap(ret, call) if ret else call
    if ret and overrides.return_shift(e):
        result = T.shift_up(result)
    sig = julia_params(e["args"], overrides.arg_shifts(e), lead=f"obj::{e['record']}")
    jl.append(f"{julia_name(name)}({sig}) = {result}")
    # `size` here is a scalar element count, so mirror it as `length` -- iterable
    # types (blocks, ProductState) need it for collect / comprehensions, since a
    # type with `iterate` defaults to IteratorSize == HasLength().
    if name == "size" and not e["args"]:
        jl.append(f"Base.length(obj::{e['record']}) = size(obj)")


# Real/complex pairs (foo / fooC): merged into one Julia function that picks the
# variant at runtime via isreal. Filled in main(); C++ keeps both registrations.
COMPLEX_BASE = {}       # (base_name, argtypes) -> C-variant entry
COMPLEX_C_SIGS = set()  # (C_name, argtypes) whose Julia is provided by the base
MERGED_NAMES = set()    # names owned by a merged pair -> same-named methods are dropped
ISREAL_TYPES = set()    # type keys with an isreal method (State, OpSum, blocks, ...)


def emit_merged(base, cvar, jl):
    args = base["args"]
    callargs = julia_call_args(args, overrides.arg_shifts(base))
    checks = [f"isreal({a['name']})" for a in args if T.decay(a["type"]) in ISREAL_TYPES]
    cond = " && ".join(checks) or "true"
    real = T.wrap(base.get("returns"), f"fun_{base['name']}({callargs})")
    cplx = T.wrap(cvar.get("returns"), f"fun_{cvar['name']}({callargs})")
    jl += [
        f"function {julia_name(base['name'])}({julia_params(args, overrides.arg_shifts(base))})",
        f"    {cond} ? ({real}) : ({cplx})",
        f"end",
    ]


def bind_function(e, cpp, jl, method_names):
    name, ret = e["name"], e.get("returns")
    names = ", ".join(a["name"] for a in e["args"])
    for params in cpp_variants(e["args"]):
        cpp.append(register_method(f"fun_{name}", params, f"{name}({names})", ret is None))
    sig = (name, tuple(a["type"] for a in e["args"]))
    if sig in COMPLEX_C_SIGS:                   # C-variant: C++ only, Julia via base
        return
    if sig in COMPLEX_BASE:                     # real base: emit merged real/complex fn
        emit_merged(e, COMPLEX_BASE[sig], jl)
        return
    if name in method_names:                    # keep C++ reg; avoid Julia clash
        return
    call = f"fun_{name}({julia_call_args(e['args'], overrides.arg_shifts(e))})"
    result = T.wrap(ret, call) if ret else call
    if ret and overrides.return_shift(e):
        result = T.shift_up(result)
    # keep the Julia-owned arrays alive while C++ views them (see csr_view)
    sparse = [a["name"] for a in e["args"] if T.category(a["type"]) == "sparse"]
    if sparse:
        result = f"GC.@preserve {' '.join(sparse)} ({result})"
    shifts = overrides.arg_shifts(e)
    jl.append(f"{julia_name(name)}({julia_params(e['args'], shifts)}) = {result}")
    if name == "to_string":
        jl.append(f"Base.show(io::IO, {julia_params(e['args'], shifts)}) = print(io, {call})")


def bind_binary_operator(e, cpp, jl):
    sym = BINARY_OPS[e["name"]]
    abi = f"op_{OP_ABI[sym]}"
    member = e.get("record") is not None
    ret, args = e.get("returns"), e["args"]
    names = ", ".join(a["name"] for a in args)
    for params in cpp_variants(args):
        full = ([cpp_self(e)] + params) if member else params
        expr = f"self {sym} {names}" if member else names.replace(", ", f" {sym} ")
        cpp.append(register_method(abi, full, expr, False))
    if member:
        b = args[0]
        call = f"{abi}(self.cxx_object, {T.unwrap(b['type'], b['name'])})"
        jl.append(f"{OP_JULIA[sym]}(self::{e['record']}, {b['name']}::{T.julia_type(b['type'])}) "
                  f"= {T.wrap(ret, call)}")
    else:
        call = f"{abi}({julia_call_args(args, set())})"
        jl.append(f"{OP_JULIA[sym]}({julia_params(args)}) = {T.wrap(ret, call)}")


def bind_index_operator(e, cpp, jl):
    ret, idx = e.get("returns"), e["args"][0]
    shifts = overrides.arg_shifts(e)
    idx_value = T.shift_down(idx["name"]) if idx["name"] in shifts else idx["name"]
    idx_jl = T.julia_type(idx["type"])
    scalar = "Scalar" in ret                    # coupling accessor (op[key] = number)
    if e.get("const"):                          # read: getindex
        if scalar:
            return   # reading a Scalar coupling needs a Scalar->number path; skip
        for params in cpp_variants(e["args"]):
            cpp.append(register_method("op_getindex", [f"{e['record']} const& self"] + params,
                                       f"self[{idx['name']}]", False))
        result = T.wrap(ret, f"op_getindex(self.cxx_object, {T.unwrap(idx['type'], idx_value)})")
        if overrides.return_shift(e):
            result = T.shift_up(result)
        jl.append(f"Base.getindex(self::{e['record']}, {idx['name']}::{idx_jl}) = {result}")
    else:                                       # write: setindex!
        val_types = ["double", "complex"] if scalar else [ret]   # Scalar accepts real/complex
        for vt in val_types:
            for params in cpp_variants(e["args"]):
                cpp.append(register_method("op_setindex", [f"{e['record']}& self"] + params +
                                           [T.cpp_param(vt, "val")], f"self[{idx['name']}] = val", True))
            jl.append(f"Base.setindex!(self::{e['record']}, val::{T.julia_type(vt)}, "
                      f"{idx['name']}::{idx_jl}) = op_setindex(self.cxx_object, "
                      f"{T.unwrap(idx['type'], idx_value)}, {T.unwrap(vt, 'val')})")


def bind_fields(e, cpp, jl):
    # Only the C++ accessor is registered; the Julia side materialises these as
    # native struct fields (see emit_types.emit), read by its `convert`.
    for f in e["fields"]:
        if not T.supported(f["type"], as_return=True):
            continue
        cpp.append(f'mod.method("mth_{f["name"]}", []({e["name"]} const& self) '
                   f'{{ return self.{f["name"]}; }});')


# --- filtering + driver -----------------------------------------------------
def wrappable(e):
    """(ok, reason). reason is a coverage bucket when not ok."""
    if e.get("record") in overrides.IGNORED_CLASSES:
        return False, "ignored class"
    reason = overrides.skip_reason(e)
    if reason:
        return False, reason
    for a in e.get("args", []):
        if not T.supported(a["type"]):
            return False, f"unsupported arg: {a['type']}"     # raw type shows the & / *
    ret = e.get("returns")
    # op[key] returning Scalar& is the coupling accessor: handled specially
    # (setindex! takes a real/complex value), so it bypasses the return gate.
    scalar_index = e["name"] == "operator[]" and ret and "Scalar" in ret
    if ret and not scalar_index and not T.supported(ret, as_return=True):
        return False, f"unsupported return: {ret}"
    name, nargs = e["name"], len(e.get("args", []))
    if name in BINARY_OPS:                       # member: self OP rhs; free: a OP b
        need = 1 if e.get("record") else 2
        if nargs != need:
            return False, f"operator arity: {name}"
    elif name == "operator[]":
        if nargs != 1 or not e.get("record"):
            return False, f"operator arity: {name}"
    elif name.startswith("operator"):
        return False, f"operator not wrapped: {name}"
    return True, None


def bind(e, cpp, jl, method_names):
    name = e["name"]
    if e["kind"] == "constructor":
        bind_constructor(e, cpp, jl)
    elif name == "operator[]":
        bind_index_operator(e, cpp, jl)
    elif name in BINARY_OPS:
        bind_binary_operator(e, cpp, jl)
    elif e["kind"] == "method":
        bind_method(e, cpp, jl)
    elif e["kind"] == "function":
        bind_function(e, cpp, jl, method_names)


def bind_cxx_iteration(typ, element, cpp, jl):
    """Iterate via the C++ iterator protocol; the iterator is mutated in place
    and the block is re-consulted for end() each step (so no iterator copy)."""
    it = f"{typ}::iterator_t"
    cpp += [
        register_method("iterate_begin", [f"{typ} const& b"], "b.begin()", False),
        register_method("iterate_end", [f"{typ} const& b"], "b.end()", False),
        register_method("iterate_incr", [f"{it}& it"], "++it", True),
        register_method("iterate_deref", [f"{it} const& it"], "*it", False),
        register_method("iterate_eq", [f"{it} const& a", f"{it} const& b"], "a == b", False),
    ]
    jl += [
        f"function Base.iterate(block::{typ})",
        f"    size(block) > 0 || return nothing",
        f"    it = iterate_begin(block.cxx_object)",
        f"    return ({element}(iterate_deref(it)), it)",
        f"end",
        f"function Base.iterate(block::{typ}, it)",
        f"    iterate_incr(it)",
        f"    iterate_eq(it, iterate_end(block.cxx_object)) ? nothing : "
        f"({element}(iterate_deref(it)), it)",
        f"end",
    ]


def bind_indexed_iteration(typ, jl):
    """Iterate by size() + getindex -- pure Julia, reuses wrapped methods.
    (Base.length comes from the `size` method binding, see bind_method.)"""
    jl += [
        f"function Base.iterate(obj::{typ}, i::Int64 = 1)",
        f"    i > size(obj) ? nothing : (obj[i], i + 1)",
        f"end",
    ]


def main():
    entries = apidb.load()
    T.register_classes(e["name"] for e in entries if e["kind"] in ("class", "struct")
                       and e["name"] not in overrides.IGNORED_CLASSES)
    T.register_result_structs(e["name"] for e in entries if e.get("fields"))
    type_module = {e["name"]: e["module"] for e in entries if e["kind"] in ("class", "struct")}
    for problem in overrides.validate(entries):
        print("OVERRIDE WARNING:", problem)

    # detect real/complex function pairs (foo / fooC with identical arg types)
    ISREAL_TYPES.update(e["record"] for e in entries
                        if e["kind"] == "method" and e["name"] == "isreal")
    if set(T.BLOCKS) <= ISREAL_TYPES:
        ISREAL_TYPES.add("Block")
    funcsig = {(e["name"], tuple(a["type"] for a in e["args"])): e for e in entries
               if e["kind"] == "function" and not e["name"].startswith("operator")}
    for (nm, at), e in funcsig.items():
        if nm.endswith("C") and (nm[:-1], at) in funcsig:
            COMPLEX_BASE[(nm[:-1], at)] = e
            COMPLEX_C_SIGS.add((nm, at))
            MERGED_NAMES.update((nm[:-1], nm))   # drop base+C methods of the pair

    by_module = defaultdict(list)
    for e in entries:
        by_module[e["module"]].append(e)

    skipped = defaultdict(list)
    n_wrapped = 0
    for module, ents in by_module.items():
        cpp, jl = [], []
        method_names = {e["name"] for e in ents if e["kind"] == "method" and wrappable(e)[0]}
        for e in ents:
            if e["kind"] in ("class", "struct"):
                if e["kind"] == "struct" and e["name"] not in overrides.IGNORED_CLASSES:
                    bind_fields(e, cpp, jl)
                continue
            ok, reason = wrappable(e)
            if not ok:
                skipped[reason].append(f"{e.get('record') or ''}::{e['name']}")
                continue
            bind(e, cpp, jl, method_names)
            n_wrapped += 1
        for typ, element in overrides.ITERABLES.items():
            if type_module.get(typ) == module:
                bind_cxx_iteration(typ, element, cpp, jl)
        for typ in overrides.INDEXED_ITERABLES:
            if type_module.get(typ) == module:
                bind_indexed_iteration(typ, jl)
        if cpp:
            emitutil.write_cpp(module, cpp)
        if jl:
            emitutil.write_jl(module, jl)

    report = [f"WRAPPED: {n_wrapped} symbols",
              f"SKIPPED: {sum(len(v) for v in skipped.values())}", ""]
    for reason, syms in sorted(skipped.items(), key=lambda kv: -len(kv[1])):
        report.append(f"  {len(syms):3d}  {reason}")
        for s in sorted(set(syms)):
            report.append(f"          {s}")
    (emitutil.OUT).mkdir(exist_ok=True)
    (emitutil.OUT / "coverage.txt").write_text("\n".join(report) + "\n")
    print("\n".join(report[:2]) + "\n-> out/coverage.txt")


if __name__ == "__main__":
    main()
