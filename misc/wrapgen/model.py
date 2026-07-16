#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
#
# libclang AST -> structured JSON model (api_model.json).
#
# Flat capture: for every XDIAG_API-marked class/struct we record its public
# members directly; free functions are captured on their own. Templates are
# ignored. Types are classified on the WRITTEN spelling (int64_t stays
# int64_t), matching the vocabulary used throughout the xdiag headers; the
# canonical type is consulted only to decide which class a wrapped type is.
#
# Speed: we parse ONE translation unit that includes every marked header (not
# one parse per header), skip function bodies, and only walk the xdiag
# namespace.

import json
import re
import shlex
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path

import clang.cindex as cindex
import apidb

PROJECT_ROOT = Path("/Users/awietek/Research/Software/xdiag")

# Written spelling -> table-key string. Anything not here (and not Block or a
# wrapped class) becomes the normalized spelling and is later reported as
# unsupported. Categorising a key lives in types.py, not here.
VOCABULARY = {
    "int64_t": "int64_t", "double": "double", "complex": "complex", "bool": "bool",
    "std::string": "std::string",
    "arma::vec": "arma::vec", "arma::cx_vec": "arma::cx_vec",
    "arma::mat": "arma::mat", "arma::cx_mat": "arma::cx_mat",
    "arma::Col<int64_t>": "arma::ivec", "arma::Mat<int64_t>": "arma::imat",
}


def normalize(spelling):
    s = spelling.replace("&&", "").replace("&", "")
    s = re.sub(r"\b(const|volatile)\b", "", s)
    return re.sub(r"\s+", " ", s).strip().replace("> >", ">>")


def base_key(t):
    """Decayed, vocab-normalized name of a value type (no const/ref/*)."""
    name = normalize(t.spelling).replace("xdiag::", "")
    if name.startswith(("std::tuple<", "std::pair<")):     # keep element types
        n = t.get_num_template_arguments()
        if n > 0:
            elems = [base_key(t.get_template_argument_type(i)) for i in range(n)]
            return "tuple<" + ", ".join(elems) + ">"
    if name in VOCABULARY:
        return VOCABULARY[name]
    if name == "Block":
        return "Block"
    for s in ("CSRMatrix", "COOMatrix", "CSCMatrix"):   # sparse family -> base name
        if name == s or name.startswith(s + "<"):
            return s
    decl = t.get_canonical().get_declaration()      # resolve which class this is
    if decl is not None and decl.spelling and decl.kind in (
            cindex.CursorKind.CLASS_DECL, cindex.CursorKind.STRUCT_DECL):
        return decl.spelling
    return name


def classify(t):
    """A clang Type -> a type string that keeps the passing qualifier, e.g.
    'int64_t', 'Op const&', 'State&', 'OpSum&&', 'double*'. types.py decays
    this to the key and uses the qualifier to decide what can be wrapped."""
    if t.kind == cindex.TypeKind.VOID:
        return "void"
    if t.kind == cindex.TypeKind.POINTER:
        return base_key(t.get_pointee()) + "*"
    if t.kind in (cindex.TypeKind.LVALUEREFERENCE, cindex.TypeKind.RVALUEREFERENCE):
        pointee = t.get_pointee()
        base = base_key(pointee)
        if t.kind == cindex.TypeKind.RVALUEREFERENCE:
            return base + "&&"
        return base + " const&" if pointee.is_const_qualified() else base + "&"
    return base_key(t)


# --- build config, file helpers ---------------------------------------------
def compile_flags():
    """The project's real compile flags, minus the compiler, the source file
    and -o <output>. Keeping the rest (crucially -isysroot / -stdlib on macOS)
    is what lets int64_t / complex / arma::* resolve instead of silently
    error-recovering to `int`."""
    entry = json.load(open(PROJECT_ROOT / "build" / "compile_commands.json"))[0]
    tokens = shlex.split(entry["command"])
    flags, skip = [], False
    for tok in tokens[1:]:                       # tokens[0] is the compiler
        if skip:                                 # drop the argument of -o
            skip = False
        elif tok == "-o":
            skip = True
        elif tok != entry["file"]:
            flags.append(tok)
    return flags


@lru_cache(maxsize=1)
def resource_dir():
    cc = shutil.which("clang") or shutil.which("clang++")
    return subprocess.run([cc, "-print-resource-dir"], capture_output=True,
                          text=True, check=True).stdout.strip()


def in_project(path):
    try:
        Path(path).resolve().relative_to(PROJECT_ROOT / "xdiag")
        return True
    except (ValueError, TypeError):
        return False


@lru_cache(maxsize=None)
def read_file(path):
    return Path(path).read_text()


def source_text(node):
    s = read_file(str(node.location.file))
    return re.sub(r"\s+", " ", s[node.extent.start.offset:node.extent.end.offset]).strip()


def module_of(header):
    parts = header.split("/")
    return parts[1] if len(parts) >= 3 else "misc"


# --- cursor -> entry (compact string form for the TOML model) ---------------
def make_arg(a):
    text = source_text(a)
    s = f"{classify(a.type)} {a.spelling}"
    if "=" in text:
        s += " = " + text.split("=", 1)[1].strip()
    return s


def is_deleted(cursor):
    try:
        return cursor.is_deleted_method()
    except Exception:
        return source_text(cursor).rstrip().endswith("= delete")


def params(cursor):
    # get_arguments() is empty for FUNCTION_TEMPLATE; fall back to PARM_DECLs
    args = list(cursor.get_arguments())
    return args or [c for c in cursor.get_children()
                    if c.kind == cindex.CursorKind.PARM_DECL]


def make_callable(cursor, kind, record, module, header):
    entry = {"kind": kind, "name": cursor.spelling, "module": module,
             "header": header, "line": cursor.location.line}
    if record:
        entry["record"] = record
    if is_deleted(cursor):
        entry["deleted"] = True
    if kind == "method" and cursor.is_const_method():
        entry["const"] = True
    args = [make_arg(a) for a in params(cursor)]
    if args:
        entry["args"] = args
    if kind in ("method", "function"):
        ret = classify(cursor.result_type)
        if ret != "void":
            entry["returns"] = ret
    return entry


def make_record(cls, module, header):
    is_struct = cls.kind == cindex.CursorKind.STRUCT_DECL
    entry = {"kind": "struct" if is_struct else "class", "name": cls.spelling,
             "module": module, "header": header, "line": cls.location.line}
    if is_struct:
        fields = [f"{classify(f.type)} {f.spelling}" for f in cls.get_children()
                  if f.kind == cindex.CursorKind.FIELD_DECL
                  and f.access_specifier == cindex.AccessSpecifier.PUBLIC]
        if fields:
            entry["fields"] = fields
    yield entry
    for member in cls.get_children():
        if member.access_specifier != cindex.AccessSpecifier.PUBLIC:
            continue
        if member.kind == cindex.CursorKind.CONSTRUCTOR:
            yield make_callable(member, "constructor", cls.spelling, module, header)
        elif member.kind == cindex.CursorKind.CXX_METHOD:      # includes operators
            yield make_callable(member, "method", cls.spelling, module, header)


def walk(node, marked_lines, entries, inside_xdiag=False):
    for c in node.get_children():
        if c.kind == cindex.CursorKind.NAMESPACE:
            if inside_xdiag or c.spelling == "xdiag":
                walk(c, marked_lines, entries, inside_xdiag=True)
            continue
        if not inside_xdiag or c.location.file is None:
            continue
        path = str(Path(c.location.file.name).resolve())
        if not in_project(path):
            continue
        # range hit: XDIAG_API may sit on a later line than the cursor start
        # (e.g. `template<...>` line 1, `XDIAG_API ret f(...)` line 2)
        lines = marked_lines.get(path, ())
        if not any(c.extent.start.line <= l <= c.extent.end.line for l in lines):
            continue
        header = str(Path(path).relative_to(PROJECT_ROOT))
        if c.kind in (cindex.CursorKind.CLASS_DECL, cindex.CursorKind.STRUCT_DECL):
            entries.extend(make_record(c, module_of(header), header))
        elif c.kind == cindex.CursorKind.FUNCTION_DECL:
            entries.append(make_callable(c, "function", None, module_of(header), header))
        elif c.kind == cindex.CursorKind.FUNCTION_TEMPLATE:
            # only sparse-family templates (the CSRMatrix algorithms); no other
            # generic templates are wrapped
            e = make_callable(c, "function", None, module_of(header), header)
            if any("CSRMatrix" in a for a in e.get("args", [])):
                entries.append(e)


def main():
    headers = [h for h in (PROJECT_ROOT / "xdiag").glob("**/*.hpp")]
    marked_lines = {}
    api_headers = []
    for h in headers:
        lines = {i for i, l in enumerate(h.read_text().splitlines(), 1) if "XDIAG_API" in l}
        if lines:
            marked_lines[str(h.resolve())] = lines
            api_headers.append(h)

    # Single translation unit including every API header -> one parse.
    umbrella = "\n".join(f'#include "{h.resolve()}"' for h in api_headers)
    args = compile_flags() + [f"-resource-dir={resource_dir()}"]
    tu = cindex.Index.create().parse(
        "xdiag_api.cpp", args=args,
        unsaved_files=[("xdiag_api.cpp", umbrella)],
        options=cindex.TranslationUnit.PARSE_SKIP_FUNCTION_BODIES)

    errors = [d for d in tu.diagnostics if d.severity >= cindex.Diagnostic.Error]
    if errors:
        print(f"WARNING: {len(errors)} parse error(s) -- types may misresolve. First few:")
        for d in errors[:3]:
            print("  ", d.spelling)

    entries = []
    walk(tu.cursor, marked_lines, entries)

    seen, unique = set(), []
    for e in entries:
        k = (e["kind"], e.get("record"), e["name"], e["header"], e["line"])
        if k not in seen:
            seen.add(k)
            unique.append(e)
    n_modules = apidb.write(unique)
    print(f"parsed {len(api_headers)} API headers, wrote {len(unique)} entries "
          f"across {n_modules} modules -> api/")


if __name__ == "__main__":
    main()
