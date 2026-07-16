#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
#
# Auto-generate a Markdown reference of the public xdiag C++ API.
#
# The public API is exactly the set of declarations annotated with the
# XDIAG_API macro.  We parse the headers with libclang, redefining XDIAG_API
# to an `annotate("xdiag_api")` attribute so the parser tags precisely those
# declarations.  We then walk the AST and emit signatures grouped by
# namespace, class and header file.
#
# Usage:
#   misc/.venv/bin/python misc/generate_api_md.py [-o OUTPUT.md]
#
# Run from anywhere; paths are resolved relative to the repository root.

import argparse
import os
import re
import shutil
import subprocess
import sys

import clang.cindex as cindex

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Marker injected in place of XDIAG_API so libclang tags the exported decls.
ANNOTATION = "xdiag_api"

# Kinds whose signature we render from raw source (callables, variables).
CALLABLE_KINDS = {
    cindex.CursorKind.FUNCTION_DECL,
    cindex.CursorKind.CXX_METHOD,
    cindex.CursorKind.CONSTRUCTOR,
    cindex.CursorKind.DESTRUCTOR,
    cindex.CursorKind.FUNCTION_TEMPLATE,
    cindex.CursorKind.CONVERSION_FUNCTION,
    cindex.CursorKind.VAR_DECL,
    cindex.CursorKind.FIELD_DECL,
}

RECORD_KINDS = {
    cindex.CursorKind.CLASS_DECL,
    cindex.CursorKind.STRUCT_DECL,
    cindex.CursorKind.CLASS_TEMPLATE,
}

_file_cache = {}


def read_file(path):
    if path not in _file_cache:
        with open(path, "r", errors="replace") as f:
            _file_cache[path] = f.read()
    return _file_cache[path]


def in_project(path):
    if not path:
        return False
    path = os.path.abspath(path)
    return (
        path.startswith(REPO_ROOT)
        and (os.sep + "xdiag" + os.sep) in path
        and (os.sep + "extern" + os.sep) not in path
        and (os.sep + "old" + os.sep) not in path
    )


def has_api_annotation(cursor):
    for child in cursor.get_children():
        if child.kind == cindex.CursorKind.ANNOTATE_ATTR and child.spelling == ANNOTATION:
            return True
    return False


def namespace_path(cursor):
    """Fully-qualified namespace of a cursor (excluding its own name)."""
    parts = []
    parent = cursor.semantic_parent
    while parent is not None and parent.kind != cindex.CursorKind.TRANSLATION_UNIT:
        if parent.kind == cindex.CursorKind.NAMESPACE and parent.spelling:
            parts.append(parent.spelling)
        parent = parent.semantic_parent
    return "::".join(reversed(parts))


def enclosing_record(cursor):
    parent = cursor.semantic_parent
    if parent is not None and parent.kind in RECORD_KINDS:
        return parent
    return None


def rel_header(cursor):
    f = cursor.location.file
    if f is None:
        return None
    return os.path.relpath(os.path.abspath(f.name), REPO_ROOT)


def raw_signature(cursor):
    """Faithful signature taken from the source text, body stripped."""
    ext = cursor.extent
    if ext.start.file is None:
        return cursor.displayname or cursor.spelling
    src = read_file(ext.start.file.name)
    text = src[ext.start.offset : ext.end.offset]
    # Cut bodies/trailing junk. Only strip `=` initializers for variables —
    # for functions `=` belongs to operator== / default arguments / = default.
    seps = ["{", ";"]
    if cursor.kind in (cindex.CursorKind.VAR_DECL, cindex.CursorKind.FIELD_DECL):
        seps.append("=")
    cuts = [i for i in (text.find(s) for s in seps) if i != -1]
    if cuts:
        text = text[: min(cuts)]
    text = re.sub(r"\bXDIAG_API\b\s*", "", text)
    text = re.sub(r"\s+", " ", text).strip().rstrip(",").strip()
    return text


def record_signature(cursor):
    kw = {
        cindex.CursorKind.CLASS_DECL: "class",
        cindex.CursorKind.STRUCT_DECL: "struct",
        cindex.CursorKind.CLASS_TEMPLATE: "class template",
    }[cursor.kind]
    return f"{kw} {cursor.spelling}"


def signature(cursor):
    if cursor.kind in RECORD_KINDS:
        return record_signature(cursor)
    if cursor.kind == cindex.CursorKind.ENUM_DECL:
        return f"enum {cursor.spelling}"
    if cursor.kind in (
        cindex.CursorKind.TYPEDEF_DECL,
        cindex.CursorKind.TYPE_ALIAS_DECL,
        cindex.CursorKind.TYPE_ALIAS_TEMPLATE_DECL,
    ):
        return raw_signature(cursor)
    return raw_signature(cursor)


def collect(cursor, out):
    """Recursively gather annotated decls, pruning out-of-project subtrees."""
    for child in cursor.get_children():
        floc = child.location.file
        # Prune subtrees living outside the project (arma, std, ...), but
        # always descend into namespaces (they may reopen in our headers).
        if floc is not None and not in_project(floc.name):
            if child.kind != cindex.CursorKind.NAMESPACE:
                continue
        if has_api_annotation(child) and in_project(floc.name if floc else None):
            out.append(child)
        collect(child, out)


def detect_system_includes():
    """Ask the system C++ compiler for its header search paths so libclang
    (which ships no toolchain of its own) can find <optional>, <vector>, ..."""
    args = []
    sdk = None
    if shutil.which("xcrun"):
        try:
            sdk = subprocess.check_output(
                ["xcrun", "--show-sdk-path"], text=True, stderr=subprocess.DEVNULL
            ).strip()
        except subprocess.CalledProcessError:
            sdk = None
    if sdk:
        args += ["-isysroot", sdk]

    compiler = (
        shutil.which("/opt/homebrew/opt/llvm/bin/clang++")
        or shutil.which("clang++")
        or shutil.which("c++")
    )
    if compiler:
        try:
            out = subprocess.run(
                [compiler, "-E", "-x", "c++", "-v", os.devnull],
                capture_output=True,
                text=True,
            ).stderr
        except OSError:
            out = ""
        capture = False
        for line in out.splitlines():
            if "#include <...> search starts here:" in line:
                capture = True
                continue
            if "End of search list." in line:
                break
            if capture:
                path = line.strip()
                if path.endswith("(framework directory)"):
                    args += ["-iframework", path.split(" (")[0]]
                elif path:
                    args += ["-isystem", path]
    return args


def parse(extra_args):
    index = cindex.Index.create()
    # Synthetic umbrella TU including every header that mentions XDIAG_API.
    headers = []
    for dirpath, _dirs, files in os.walk(os.path.join(REPO_ROOT, "xdiag")):
        if os.sep + "old" + os.sep in dirpath + os.sep:
            continue
        for fn in files:
            if not fn.endswith(".hpp") or fn == "all.hpp":
                continue
            full = os.path.join(dirpath, fn)
            if "XDIAG_API" in read_file(full):
                headers.append(os.path.relpath(full, REPO_ROOT))
    headers.sort()

    umbrella = "\n".join(f'#include "{h}"' for h in headers) + "\n"
    # Override xdiag_api.hpp so XDIAG_API is always our annotate attribute.
    # (The real header unconditionally re-#defines XDIAG_API, which would
    # otherwise clobber a command-line -DXDIAG_API.)
    api_hpp = os.path.join(REPO_ROOT, "xdiag", "utils", "xdiag_api.hpp")
    api_stub = (
        "#pragma once\n"
        f'#define XDIAG_API __attribute__((annotate("{ANNOTATION}")))\n'
    )
    unsaved = [("__xdiag_api__.cpp", umbrella), (api_hpp, api_stub)]

    args = [
        "-x",
        "c++",
        "-std=c++17",
        f"-I{REPO_ROOT}",
        "-I/opt/homebrew/include",  # hdf5.h
        "-DXDIAG_STATIC_DEFINE",
        "-DXDIAG_USE_OPENMP=0",
        "-ferror-limit=0",
    ] + detect_system_includes() + extra_args

    tu = index.parse(
        "__xdiag_api__.cpp",
        args=args,
        unsaved_files=unsaved,
        options=cindex.TranslationUnit.PARSE_SKIP_FUNCTION_BODIES
        | cindex.TranslationUnit.PARSE_INCOMPLETE,
    )
    # Surface fatal diagnostics but keep going — the annotate attribute is
    # attached regardless of unresolved types.
    fatals = [d for d in tu.diagnostics if d.severity >= cindex.Diagnostic.Error]
    return tu, fatals, len(headers)


def render(cursors):
    # Group: namespace -> {"free": [...], "records": {record_name: [...]}}
    namespaces = {}
    for c in cursors:
        ns = namespace_path(c)
        rec = enclosing_record(c)
        grp = namespaces.setdefault(ns, {"free": [], "records": {}})
        entry = {
            "sig": signature(c),
            "file": rel_header(c),
            "line": c.location.line,
            "kind": c.kind,
        }
        if rec is not None:
            # Namespace of the record itself, not the member.
            grp["records"].setdefault(rec.spelling, []).append(entry)
        elif c.kind in RECORD_KINDS:
            grp["records"].setdefault(c.spelling, []).append({**entry, "is_decl": True})
        else:
            grp["free"].append(entry)

    lines = []
    lines.append("# xdiag C++ API reference")
    lines.append("")
    lines.append(
        "_Auto-generated by `misc/generate_api_md.py` from declarations marked "
        "`XDIAG_API`. Do not edit by hand._"
    )
    lines.append("")

    for ns in sorted(namespaces, key=lambda s: (s == "", s)):
        grp = namespaces[ns]
        title = f"namespace `{ns}`" if ns else "global namespace"
        lines.append(f"## {title}")
        lines.append("")

        # Free functions / variables grouped by header.
        by_file = {}
        for e in grp["free"]:
            by_file.setdefault(e["file"], []).append(e)
        for fpath in sorted(by_file):
            lines.append(f"### `{fpath}`")
            lines.append("")
            for e in sorted(by_file[fpath], key=lambda e: e["line"]):
                lines.append(f"- `{e['sig']}`")
            lines.append("")

        # Records (classes/structs) and their members.
        for rname in sorted(grp["records"]):
            entries = grp["records"][rname]
            decl = next((e for e in entries if e.get("is_decl")), None)
            members = [e for e in entries if not e.get("is_decl")]
            files = sorted({e["file"] for e in entries})
            qual = f"{ns}::{rname}" if ns else rname
            kw = decl["sig"].split()[0] if decl else "class"
            lines.append(f"### {kw} `{qual}`")
            lines.append("")
            lines.append(", ".join(f"`{f}`" for f in files))
            lines.append("")
            for e in sorted(members, key=lambda e: e["line"]):
                lines.append(f"- `{e['sig']}`")
            lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "-o",
        "--output",
        default="-",
        help="output markdown file (default: stdout)",
    )
    ap.add_argument(
        "clang_args",
        nargs="*",
        help="extra args forwarded to clang (e.g. -I paths)",
    )
    args = ap.parse_args()

    tu, fatals, n_headers = parse(args.clang_args)
    cursors = []
    collect(tu.cursor, cursors)

    md = render(cursors)
    if args.output == "-":
        sys.stdout.write(md)
    else:
        out = (
            args.output
            if os.path.isabs(args.output)
            else os.path.join(os.getcwd(), args.output)
        )
        with open(out, "w") as f:
            f.write(md)
        print(
            f"Wrote {out}: {len(cursors)} API symbols from {n_headers} headers.",
            file=sys.stderr,
        )

    if fatals:
        print(
            f"\nNote: {len(fatals)} parse error(s) (unresolved includes/types). "
            "API extraction is unaffected unless symbols are missing. "
            "Pass include paths as extra args if needed, e.g.:\n"
            "  misc/generate_api_md.py -- -I/path/to/hdf5/include",
            file=sys.stderr,
        )
        for d in fatals[:5]:
            print(f"  {d.spelling}", file=sys.stderr)


if __name__ == "__main__":
    main()
