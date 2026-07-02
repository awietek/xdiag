# SPDX-License-Identifier: Apache-2.0
#
# Structured model of the public xdiag C++ API for Julia wrapper generation.
#
# Reuses the libclang machinery in ../generate_api_md.py to find every
# declaration marked XDIAG_API, then builds a structured, JSON-serialisable
# model in which every argument and return type is *classified* into the
# categories the wrapper generator cares about (primitive / string / arma /
# std::vector / wrapped-type / raw-pointer / ...).
#
# Type classification is done on the *canonical* clang type, never on the
# display spelling, so typedefs (arma::vec, xdiag::complex, per-class
# iterator_t, ...) are resolved consistently.
#
# Usage:
#   misc/.venv/bin/python misc/wrapgen/model.py --json misc/wrapgen/api_model.json
#   misc/.venv/bin/python misc/wrapgen/model.py --report -   # coverage report

import argparse
import json
import os
import re
import sys

MISC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(MISC_DIR)
sys.path.insert(0, MISC_DIR)

import clang.cindex as cindex  # noqa: E402
import generate_api_md as api  # noqa: E402  (reuse parse/collect/helpers)

# ---------------------------------------------------------------------------
# Type classification
# ---------------------------------------------------------------------------

_ARMA_CANON = {
    "arma::Col<double>": "vec",
    "arma::Col<std::complex<double>>": "cx_vec",
    "arma::Mat<double>": "mat",
    "arma::Mat<std::complex<double>>": "cx_mat",
    "arma::Col<long long>": "ivec",
    "arma::Mat<long long>": "imat",
}

_PRIMITIVE_CANON = {
    "bool": "bool",
    "char": "char",
    "int": "int32",
    "long": "int64",
    "long long": "int64",
    "unsigned long": "uint64",
    "unsigned long long": "uint64",
    "unsigned int": "uint32",
    "short": "int16",
    "double": "double",
    "float": "float",
    "std::complex<double>": "complex",
}

_STRING_CANON = {
    "std::string",
    "std::basic_string<char>",
    "const char *",
    "const char*",
    "char *",
}


def _norm(spelling):
    """Normalise a canonical type spelling: drop const/volatile/&, collapse ws."""
    s = spelling
    s = s.replace("&&", "").replace("&", "")
    s = re.sub(r"\bconst\b", "", s)
    s = re.sub(r"\bvolatile\b", "", s)
    s = re.sub(r"\s+", " ", s).strip()
    # normalise "> >" -> ">>" for map lookups
    s = s.replace("> >", ">>").replace("> >", ">>")
    return s


def classify_type(clang_type, record_names):
    """Classify a clang Type into a category descriptor dict.

    record_names: set of fully-qualified record names that are wrapped types.
    """
    t = clang_type.get_canonical()

    # Pointers -> raw pointer (routed to hand-written specials)
    if t.kind == cindex.TypeKind.POINTER:
        pointee = t.get_pointee()
        inner = classify_type(pointee, record_names)
        return {"cat": "pointer", "cpp": _norm(t.spelling), "elem": inner}

    # References -> classify the pointee, keep the ref flag
    if t.kind in (cindex.TypeKind.LVALUEREFERENCE, cindex.TypeKind.RVALUEREFERENCE):
        pointee = t.get_pointee()
        inner = classify_type(pointee, record_names)
        inner = dict(inner)
        inner["ref"] = "rvalue" if t.kind == cindex.TypeKind.RVALUEREFERENCE else "lvalue"
        inner["const"] = pointee.is_const_qualified()
        return inner

    if t.kind == cindex.TypeKind.VOID:
        return {"cat": "void", "cpp": "void"}

    norm = _norm(t.spelling)
    norm_nospace = norm.replace(" ", "")

    # Unresolved template parameters (from FUNCTION_TEMPLATE / class templates):
    # these need explicit instantiation lists in the overrides.
    if "type-parameter-" in norm:
        return {"cat": "template_param", "cpp": norm}

    if norm in _PRIMITIVE_CANON:
        return {"cat": "primitive", "cpp": norm, "prim": _PRIMITIVE_CANON[norm]}

    if norm in _STRING_CANON or norm_nospace in {s.replace(" ", "") for s in _STRING_CANON}:
        return {"cat": "string", "cpp": norm}

    if norm in _ARMA_CANON:
        return {"cat": "arma", "cpp": norm, "arma": _ARMA_CANON[norm]}
    if norm_nospace in {k.replace(" ", "") for k in _ARMA_CANON}:
        key = next(k for k in _ARMA_CANON if k.replace(" ", "") == norm_nospace)
        return {"cat": "arma", "cpp": norm, "arma": _ARMA_CANON[key]}

    # std::vector<...>
    if norm.startswith("std::vector<"):
        elem = None
        try:
            targ = t.get_template_argument_type(0)
            if targ.kind != cindex.TypeKind.INVALID:
                elem = classify_type(targ, record_names)
        except Exception:
            elem = None
        return {"cat": "stdvector", "cpp": norm, "elem": elem or {"cat": "unknown"}}

    # std::variant<...> — notably the polymorphic Block (Spinhalf|Boson|...).
    if norm.startswith("std::variant<"):
        members = []
        try:
            i = 0
            while True:
                targ = t.get_template_argument_type(i)
                if targ.kind == cindex.TypeKind.INVALID:
                    break
                members.append(classify_type(targ, record_names))
                i += 1
        except Exception:
            pass
        is_block = bool(members) and all(
            mm.get("cat") == "wrapped" for mm in members)
        return {"cat": "variant", "cpp": norm, "members": members,
                "is_block": is_block}

    # Other std:: composites -> flagged for a rule / override
    if norm.startswith("std::pair<") or norm.startswith("std::tuple<") \
       or norm.startswith("std::map<") or norm.startswith("std::optional<") \
       or norm.startswith("std::__wrap_iter<") or norm.startswith("std::chrono::"):
        return {"cat": "std_composite", "cpp": norm}

    # A wrapped xdiag type? (any in-project record declaration)
    decl = t.get_declaration()
    qual = _qualified_name(decl) if decl is not None else None
    if qual and qual in record_names:
        return {"cat": "wrapped", "cpp": norm, "type": qual}

    # Enum
    if t.kind == cindex.TypeKind.ENUM:
        return {"cat": "enum", "cpp": norm, "type": qual}

    if norm.startswith("std::ostream") or norm.startswith("std::basic_ostream"):
        return {"cat": "ostream", "cpp": norm}

    return {"cat": "unknown", "cpp": norm, "type": qual}


# ---------------------------------------------------------------------------
# AST -> model
# ---------------------------------------------------------------------------

def _qualified_name(cursor):
    if cursor is None or not cursor.spelling:
        return None
    parts = [cursor.spelling]
    parent = cursor.semantic_parent
    while parent is not None and parent.kind != cindex.CursorKind.TRANSLATION_UNIT:
        if parent.spelling:
            parts.append(parent.spelling)
        parent = parent.semantic_parent
    return "::".join(reversed(parts))


_OP_RE = re.compile(r"^operator\s*(.+)$")


def _collect_all_records(cursor, out):
    """Walk the whole TU for in-project record declarations (regardless of
    annotation).  The XDIAG_API annotation propagates to *members* but is not
    reliably detectable on the class node itself, so we recognise wrapped
    types by their declaration instead.  Keyed by qualified name; a definition
    location wins over a forward declaration."""
    for child in cursor.get_children():
        floc = child.location.file
        in_proj = floc is not None and api.in_project(floc.name)
        if child.kind in api.RECORD_KINDS and in_proj and child.spelling:
            q = _qualified_name(child)
            if q and (q not in out or child.is_definition()):
                out[q] = {
                    "qualified": q,
                    "record_kind": {
                        cindex.CursorKind.CLASS_DECL: "class",
                        cindex.CursorKind.STRUCT_DECL: "struct",
                        cindex.CursorKind.CLASS_TEMPLATE: "class_template",
                    }[child.kind],
                    "header": api.rel_header(child),
                    "line": child.location.line,
                }
        # Descend everywhere in-project (and into namespaces regardless).
        if in_proj or child.kind == cindex.CursorKind.NAMESPACE:
            _collect_all_records(child, out)


def build_model(extra_args=None):
    tu, fatals, n_headers = api.parse(extra_args or [])
    cursors = []
    api.collect(tu.cursor, cursors)

    # All in-project record decls -> the set of names that classify as wrapped.
    rec_map = {}
    _collect_all_records(tu.cursor, rec_map)
    record_names = set(rec_map)

    entries = []
    for c in cursors:
        e = _entry(c, record_names)
        if e is not None:
            entries.append(e)

    # Which records are actually part of the API surface: those that enclose a
    # collected member, or are referenced as a wrapped arg/return type.
    api_record_names = set()
    for e in entries:
        if e.get("record"):
            # map bare record name back to a qualified name via rec_map
            for q in record_names:
                if q.split("::")[-1] == e["record"]:
                    api_record_names.add(q)
        for tp in [e.get("return")] + [a["type"] for a in e.get("args", [])]:
            if tp and tp.get("cat") == "wrapped":
                api_record_names.add(tp["type"])
            if tp and tp.get("cat") == "variant":
                for mm in tp.get("members", []):
                    if mm.get("cat") == "wrapped":
                        api_record_names.add(mm["type"])

    # Synthesise record entries for the report / downstream emitters.
    record_entries = []
    for q in sorted(api_record_names):
        info = rec_map.get(q)
        if info is None:
            continue
        record_entries.append({
            "kind": "record",
            "name": q.split("::")[-1],
            "qualified": q,
            "namespace": "::".join(q.split("::")[:-1]),
            "record": None,
            "header": info["header"],
            "line": info["line"],
            "record_kind": info["record_kind"],
            "signature": f"{info['record_kind'].replace('_', ' ')} {q.split('::')[-1]}",
        })

    return {
        "n_headers": n_headers,
        "n_fatals": len(fatals),
        "record_names": sorted(record_names),
        "api_record_names": sorted(api_record_names),
        "entries": record_entries + entries,
    }


def _entry(c, record_names):
    header = api.rel_header(c)
    ns = api.namespace_path(c)
    rec = api.enclosing_record(c)
    rec_name = rec.spelling if rec is not None else None
    base = {
        "name": c.spelling,
        "qualified": _qualified_name(c),
        "namespace": ns,
        "record": rec_name,
        "header": header,
        "line": c.location.line,
        "signature": re.sub(r"^[A-Za-z_]*API\b\s*", "", api.signature(c)),
    }

    if c.kind in api.RECORD_KINDS:
        base["kind"] = "record"
        base["record_kind"] = {
            cindex.CursorKind.CLASS_DECL: "class",
            cindex.CursorKind.STRUCT_DECL: "struct",
            cindex.CursorKind.CLASS_TEMPLATE: "class_template",
        }[c.kind]
        return base

    if c.kind == cindex.CursorKind.ENUM_DECL:
        base["kind"] = "enum"
        return base

    if c.kind in (cindex.CursorKind.TYPEDEF_DECL,
                  cindex.CursorKind.TYPE_ALIAS_DECL,
                  cindex.CursorKind.TYPE_ALIAS_TEMPLATE_DECL):
        base["kind"] = "typedef"
        return base

    if c.kind in (cindex.CursorKind.VAR_DECL, cindex.CursorKind.FIELD_DECL):
        base["kind"] = "var"
        base["type"] = classify_type(c.type, record_names)
        return base

    # Callables
    if c.kind in api.CALLABLE_KINDS:
        args = []
        for a in _params(c):
            entry = {
                "name": a.spelling,
                "type": classify_type(a.type, record_names),
                "has_default": _has_default(a),
            }
            if entry["has_default"]:
                dv = _default_value(base["signature"], a.spelling)
                if dv is not None:
                    entry["default"] = dv
            args.append(entry)
        base["args"] = args
        base["is_template"] = c.kind == cindex.CursorKind.FUNCTION_TEMPLATE
        tparams = _template_params(c)
        if tparams:
            base["template_params"] = tparams
        base["is_static"] = c.is_static_method() if hasattr(c, "is_static_method") else False
        try:
            base["is_const"] = c.is_const_method()
        except Exception:
            base["is_const"] = False

        if c.kind == cindex.CursorKind.CONSTRUCTOR:
            base["kind"] = "constructor"
        elif c.kind == cindex.CursorKind.DESTRUCTOR:
            base["kind"] = "destructor"
        else:
            base["return"] = classify_type(c.result_type, record_names)
            m = _OP_RE.match(c.spelling)
            if m:
                base["kind"] = "operator"
                base["op"] = m.group(1).strip()
            elif rec_name is not None:
                base["kind"] = "method"
            else:
                base["kind"] = "function"
        return base

    return None


_TEMPLATE_PARAM_KINDS = (
    cindex.CursorKind.TEMPLATE_TYPE_PARAMETER,
    cindex.CursorKind.TEMPLATE_NON_TYPE_PARAMETER,
    cindex.CursorKind.TEMPLATE_TEMPLATE_PARAMETER,
)


def _params(cursor):
    """Parameter cursors. libclang's get_arguments() is empty for
    FUNCTION_TEMPLATE, so fall back to walking PARM_DECL children."""
    args = list(cursor.get_arguments())
    if args:
        return args
    return [ch for ch in cursor.get_children()
            if ch.kind == cindex.CursorKind.PARM_DECL]


def _template_params(cursor):
    return [ch.spelling for ch in cursor.get_children()
            if ch.kind in _TEMPLATE_PARAM_KINDS]


def _has_default(arg_cursor):
    """Heuristic: does this parameter have a default argument?"""
    for tok in arg_cursor.get_tokens():
        if tok.spelling == "=":
            return True
    return False


def _default_value(signature, arg_name):
    """Extract the C++ default value for arg_name from the raw signature, e.g.
    '..., std::string backend = "auto")' -> '"auto"'. Handles the simple
    literal defaults used in the API (numbers, strings, bools)."""
    if not arg_name:
        return None
    m = re.search(r"\b" + re.escape(arg_name) + r"\s*=\s*([^,)]+)", signature)
    if not m:
        return None
    return m.group(1).strip()


# ---------------------------------------------------------------------------
# Coverage report (against the current hand-written wrapper)
# ---------------------------------------------------------------------------

def current_wrapped():
    """Best-effort scan of the hand-written C++ wrapper: which types have an
    add_type<>, and which method names appear (heuristic; renames not tracked)."""
    src_dir = os.path.join(REPO_ROOT, "julia", "src")
    types = set()
    methods = set()  # bare method-name strings seen anywhere
    for dp, _d, files in os.walk(src_dir):
        for fn in files:
            if not fn.endswith(".cpp"):
                continue
            with open(os.path.join(dp, fn), errors="replace") as f:
                txt = f.read()
            for m in re.finditer(r"add_type<\s*([A-Za-z0-9_:<>, ]+?)\s*>", txt):
                types.add(m.group(1).split("<")[0].split("::")[-1].strip())
            for m in re.finditer(r'\.method\(\s*"([^"]+)"', txt):
                methods.add(m.group(1))
            for m in re.finditer(r'mod\.method\(\s*"([^"]+)"', txt):
                methods.add(m.group(1))
    return types, methods


# Types that look like internal algebra machinery -> flag for a scope decision.
INTERNAL_CANDIDATES = {
    "Scalar", "Coeff", "Monomial", "Term", "Vector", "Matrix",
}


def report(model):
    types, wrapped_methods = current_wrapped()
    lines = []
    L = lines.append

    recs = [e for e in model["entries"] if e["kind"] == "record"]
    by_name = {}
    for e in model["entries"]:
        if e["record"]:
            by_name.setdefault(e["record"], []).append(e)

    free = [e for e in model["entries"]
            if e["kind"] in ("function", "operator") and not e["record"]]

    L("# xdiag Julia wrapper — API inventory & coverage")
    L("")
    L(f"- headers scanned: {model['n_headers']}  (parse errors: {model['n_fatals']})")
    L(f"- API records (classes/structs): {len(recs)}")
    L(f"- free functions/operators: {len(free)}")
    L(f"- currently add_type<>'d in julia/src: {len(types)}")
    L("")

    # Category tally over all callable args/returns -> shows classifier health.
    cat_tally = {}
    unknown_types = {}
    for e in model["entries"]:
        tps = []
        if "return" in e:
            tps.append(e["return"])
        for a in e.get("args", []):
            tps.append(a["type"])
        for tp in tps:
            cat_tally[tp["cat"]] = cat_tally.get(tp["cat"], 0) + 1
            if tp["cat"] == "unknown":
                unknown_types[tp["cpp"]] = unknown_types.get(tp["cpp"], 0) + 1
    L("## Type-classification tally (args + returns)")
    L("")
    for cat, n in sorted(cat_tally.items(), key=lambda kv: -kv[1]):
        L(f"- {cat}: {n}")
    L("")
    if unknown_types:
        L("### Unclassified ('unknown') types — need a classifier rule or override")
        L("")
        for cpp, n in sorted(unknown_types.items(), key=lambda kv: -kv[1]):
            L(f"- `{cpp}` ×{n}")
        L("")

    # Internal-type candidates: surface full member lists for a scope decision.
    L("## Internal-type candidates (scope decision)")
    L("")
    L("_Full member lists so you can decide wrap-vs-exclude type-by-type._")
    L("")
    for e in sorted(recs, key=lambda e: e["name"]):
        if e["name"] not in INTERNAL_CANDIDATES:
            continue
        members = by_name.get(e["name"], [])
        L(f"### `{e['qualified']}` — {len(members)} members  "
          f"({'currently wrapped' if e['name'] in types else 'NOT wrapped'})")
        L("")
        for m in sorted(members, key=lambda m: m["line"]):
            L(f"- `{m['signature']}`")
        L("")

    # Distributed (MPI) — excluded from serial Julia lib.
    L("## Distributed / MPI records (excluded from serial Julia wrapper)")
    L("")
    for e in sorted(recs, key=lambda e: e["name"]):
        if "distributed" in (e["header"] or "").lower() or "Distributed" in e["name"]:
            L(f"- `{e['qualified']}`  ({e['header']})")
    L("")

    # Coverage table for the main records.
    L("## Record coverage")
    L("")
    for e in sorted(recs, key=lambda e: e["name"]):
        if "Distributed" in e["name"]:
            continue
        members = by_name.get(e["name"], [])
        status = "✓ wrapped" if e["name"] in types else "✗ MISSING"
        L(f"### `{e['qualified']}`  [{status}] — {len(members)} members")
        L("")
        for m in sorted(members, key=lambda m: m["line"]):
            mark = ""
            op = m.get("op")
            if m["name"] in wrapped_methods or (op and op in wrapped_methods):
                mark = " ✓"
            L(f"- `{m['signature']}`{mark}")
        L("")

    return "\n".join(lines).rstrip() + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", help="write JSON model to this path")
    ap.add_argument("--report", help="write coverage report to this path ('-' = stdout)")
    ap.add_argument("clang_args", nargs="*")
    args = ap.parse_args()

    model = build_model(args.clang_args)

    if args.json:
        with open(args.json, "w") as f:
            json.dump(model, f, indent=2)
        print(f"Wrote {args.json}: {len(model['entries'])} entries.", file=sys.stderr)

    if args.report:
        md = report(model)
        if args.report == "-":
            sys.stdout.write(md)
        else:
            with open(args.report, "w") as f:
                f.write(md)
            print(f"Wrote {args.report}.", file=sys.stderr)

    if not args.json and not args.report:
        print(f"Model built: {len(model['entries'])} entries, "
              f"{len(model['record_names'])} records. "
              f"Pass --json / --report to emit.", file=sys.stderr)


if __name__ == "__main__":
    main()
