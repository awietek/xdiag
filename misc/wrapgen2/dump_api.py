#!/usr/bin/env python3
import re
import sys
import json
import re
import shlex
from pathlib import Path
from functools import lru_cache
import subprocess, shutil
import clang.cindex as cindex

API_KINDS = {
    cindex.CursorKind.FUNCTION_DECL, cindex.CursorKind.CXX_METHOD,
    cindex.CursorKind.CONSTRUCTOR, cindex.CursorKind.DESTRUCTOR,
    cindex.CursorKind.FUNCTION_TEMPLATE, cindex.CursorKind.CLASS_TEMPLATE,
    cindex.CursorKind.CLASS_DECL, cindex.CursorKind.STRUCT_DECL,
    cindex.CursorKind.ENUM_DECL,
}
LEAF_KINDS = {
    cindex.CursorKind.FUNCTION_DECL, cindex.CursorKind.CXX_METHOD,
    cindex.CursorKind.CONSTRUCTOR, cindex.CursorKind.DESTRUCTOR,
    cindex.CursorKind.FUNCTION_TEMPLATE,
}

CLASS_KINDS = {cindex.CursorKind.CLASS_DECL, cindex.CursorKind.STRUCT_DECL}

MEMBER_KINDS = {
    cindex.CursorKind.CXX_METHOD,
    cindex.CursorKind.CONSTRUCTOR,
    cindex.CursorKind.DESTRUCTOR,
    cindex.CursorKind.FUNCTION_TEMPLATE,  # templated methods
}


def get_project_flags(project_root: Path) -> list[str]:
    compile_db_path = project_root / "build" / "compile_commands.json"
    entries = json.load(open(compile_db_path))
    args = shlex.split(entries[0]["command"])[1:-1]  # drop compiler binary + filename
    return [a for a in args if a.startswith(("-I", "-D", "-std", "-isystem"))]

def is_in_project(path: Path | None, project_root: Path) -> bool:
    if path is None:
        return False
    try:
        path.resolve().relative_to(project_root / Path("xdiag"))
        return True
    except ValueError:
        return False

@lru_cache(maxsize=1)
def resource_dir():
    clang_bin = shutil.which("clang") or shutil.which("clang++")
    return subprocess.run([clang_bin, "-print-resource-dir"],
                           capture_output=True, text=True, check=True).stdout.strip()


@lru_cache(maxsize=None)
def _read_file(path: str) -> str:
    return Path(path).read_text()

def get_source_text(node) -> str:
    text = _read_file(str(node.location.file))
    return text[node.extent.start.offset:node.extent.end.offset]

def clean_signature(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def find_api_lines(header: Path) -> set[int]:
    return {i for i, l in enumerate(header.read_text().splitlines(), 1) if "XDIAG_API" in l}

def print_api_by_file(results, project_root: Path):
    by_file = {}
    for kind, name, path, line, signature, parent_class in results:
        by_file.setdefault(path, []).append((line, kind, name, signature, parent_class))
        
    with open("api.txt", "w") as fl:
        for path in sorted(by_file, key=lambda p: Path(p).relative_to(project_root).parts):
            rel = Path(path).relative_to(project_root)
            ss = f"\n=== {rel} ==="
            print(ss)
            fl.write("{}\n".format(ss))
        
            for line, kind, name, signature, parent_class in sorted(set(by_file[path])):
                if parent_class == None:
                    ss = f"[{line:20}] {kind:20s} {signature}"
                else:
                    ss = f"[{parent_class:20s}] {kind:20s} {signature}"
                print(ss)
                fl.write("{}\n".format(ss))  

def collect_api(node, project_root: Path, api_lines: dict, results: list, inherited_class: str = None):
    loc_file = node.location.file
    path = Path(loc_file.name).resolve() if loc_file else None

    if path is not None and not is_in_project(path, project_root):
        return

    resolved = str(path) if path is not None else None
    is_class = node.kind in CLASS_KINDS
    line_matched = (
        resolved is not None
        and node.extent.start.line in api_lines.get(resolved, set())
    )

    new_inherited_class = inherited_class
    if is_class:
        new_inherited_class = node.spelling if line_matched else None

    # Case 1: class itself marked XDIAG_API
    if is_class and line_matched:
        signature = clean_signature(get_source_text(node))
        results.append((node.kind.name, node.spelling, resolved, node.location.line, signature, None))

    # Case 2: standalone function/template marked directly with XDIAG_API
    elif not is_class and inherited_class is None and node.kind in API_KINDS and line_matched:
        signature = clean_signature(get_source_text(node))
        results.append((node.kind.name, node.spelling, resolved, node.location.line, signature, None))
        if node.kind in LEAF_KINDS:
            return

    # Case 3: member of an XDIAG_API-marked class
    elif (
        inherited_class is not None
        and node.kind in MEMBER_KINDS
        and node.access_specifier == cindex.AccessSpecifier.PUBLIC
    ):
        signature = clean_signature(get_source_text(node))
        results.append((node.kind.name, node.spelling, resolved, node.location.line, signature, inherited_class))
        return

    for child in node.get_children():
        collect_api(child, project_root, api_lines, results, inherited_class=new_inherited_class)
            
def main():
    project_root = "/Users/awietek/Research/Software/xdiag"
    headers = list((Path(project_root) / Path("xdiag")) .glob("**/*.hpp"))
    api_lines = {str(h.resolve()): find_api_lines(h) for h in headers}

    project_flags = get_project_flags(Path(project_root))
    parse_args = project_flags + [f"-resource-dir={resource_dir()}"]
    
    index = cindex.Index.create()
    results = []
    for header in headers:
        print("Scanning {}".format(str(header.relative_to(project_root))))

        if not api_lines[str(header.resolve())]:
            continue  # skip headers with no XDIAG_API at all — free speedup

        tu = index.parse(
            str(header), args=parse_args,
            options=cindex.TranslationUnit.PARSE_SKIP_FUNCTION_BODIES,
        )
        collect_api(tu.cursor, project_root, api_lines, results)

    print_api_by_file(results, project_root)

if __name__ == "__main__":
    main()
