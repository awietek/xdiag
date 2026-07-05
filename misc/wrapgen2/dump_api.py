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


# def get_resolved_signature(node) -> str:
#     if node.kind in (cindex.CursorKind.FUNCTION_DECL, cindex.CursorKind.CXX_METHOD,
#                       cindex.CursorKind.CONSTRUCTOR, cindex.CursorKind.DESTRUCTOR):
#         return_type = node.result_type.get_canonical().spelling if node.kind != cindex.CursorKind.CONSTRUCTOR else ""
#         params = ", ".join(f"{arg.type.spelling} {arg.spelling}".strip()
#                             for arg in node.get_arguments())
#         return f"{return_type} {node.spelling}({params})".strip()
#     return node.type.spelling

def find_api_lines(header: Path) -> set[int]:
    return {i for i, l in enumerate(header.read_text().splitlines(), 1) if "XDIAG_API" in l}

def collect_api(node, project_root: str, api_lines: dict, results: list):
    loc_file = node.location.file
    path = Path(loc_file.name).resolve() if loc_file else None


    # Only prune when we actually have a location that's outside the project.
    # A None location (e.g. the TU root) should NOT be pruned — just don't
    # try to collect it as an API node, but still descend into its children.
    if path is not None and not is_in_project(path, project_root):
        return

    if path is not None and node.kind in API_KINDS:
        if node.extent.start.line in api_lines.get(str(path), set()):
            # results.append((node.kind.name, node.spelling, str(path), node.location.line))
            results.append((node.kind.name,
                            node.spelling,                            
                            str(path), node.location.line,
                            clean_signature(get_source_text(node)),))
            if node.kind in LEAF_KINDS:
                return

    for child in node.get_children():
        collect_api(child, project_root, api_lines, results)

def print_api_by_file(results, project_root: Path):
    by_file = {}
    for kind, name, path, line, signature in results:
        by_file.setdefault(path, []).append((line, kind, name, signature))

    for path in sorted(by_file, key=lambda p: Path(p).relative_to(project_root).parts):
        rel = Path(path).relative_to(project_root)
        print(f"\n=== {rel} ===")
        for line, kind, name, signature in sorted(set(by_file[path])):
            print(f"  [{line}] {kind:20s} {signature}")
        
def main():
    project_root = str(Path(sys.argv[1]).resolve())
    headers = list((Path(project_root) / Path("xdiag")) .glob("**/*.hpp"))
    api_lines = {str(h.resolve()): find_api_lines(h) for h in headers}

    project_flags = get_project_flags(Path(project_root))
    parse_args = project_flags + [f"-resource-dir={resource_dir()}"]
    
    index = cindex.Index.create()
    results = []
    for header in headers:
        print("{}".format(str(header.relative_to(project_root))))

        if not api_lines[str(header.resolve())]:
            continue  # skip headers with no XDIAG_API at all — free speedup

        tu = index.parse(
            str(header), args=parse_args,
            options=cindex.TranslationUnit.PARSE_SKIP_FUNCTION_BODIES,
        )
        collect_api(tu.cursor, project_root, api_lines, results)

    print_api_by_file(results, project_root)
        
    # def sort_key(entry):
    #     kind, name, path, line = entry
    #     rel = Path(path).relative_to(project_root)
    #     parts = rel.parts  # e.g. ('states', 'state.hpp') or ('math', 'matrix.hpp')
    #     return (parts[:-1], parts[-1], line, name)

        
    # for kind, name, path, line in sorted(set(results), key=sort_key):
    #     # print(f"{kind:20s} {name:30s} {Path(path).relative_to(project_root)}:{line}")
    #     print(f"{name:70s} {Path(path).relative_to(project_root)}:{line}")
        

if __name__ == "__main__":
    main()
