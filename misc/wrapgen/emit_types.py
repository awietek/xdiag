#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
#
# Emit the type definitions: one add_type<> per wrapped class on the C++ side
# (out/cpp/types.cpp), and the matching Julia struct + convert on the Julia
# side (out/julia/types.jl). Blocks additionally get `<: Block`. These must be
# defined before any method, which is why they live in their own files.

import apidb
import emitutil
import overrides
import typetable


def emit(records):
    cpp = [f'mod.add_type<{r["name"]}>("typ_{r["name"]}");' for r in records]
    jl = ["abstract type Block end"]
    # Plain wrapped types first (result-struct fields may reference them, e.g.
    # eigenvectors::State), result structs (records with public fields) last.
    plain = [r for r in records if not r.get("fields")]
    results = [r for r in records if r.get("fields")]
    for r in plain:
        name = r["name"]
        supertype = " <: Block" if name in typetable.BLOCKS else ""
        jl.append(f"\nstruct {name}{supertype}\n    cxx_object::typ_{name}\nend")
        jl.append(f"Base.convert(::Type{{Q}}, obj::typ_{name}) where {{Q<:{name}}} = {name}(obj)")
    for r in results:
        name = r["name"]
        fields = [f for f in r["fields"] if typetable.supported(f["type"], as_return=True)]
        decls = "\n".join(f"    {f['name']}::{typetable.julia_type(f['type'])}" for f in fields)
        # convert pulls each field off the C++ object via its mth_ accessor.
        vals = ", ".join(typetable.wrap(f["type"], f'mth_{f["name"]}(obj)') for f in fields)
        jl.append(f"\nstruct {name}\n{decls}\nend")
        jl.append(f"Base.convert(::Type{{Q}}, obj::typ_{name}) where {{Q<:{name}}} = "
                  f"{name}({vals})")
    return cpp, jl


def main():
    entries = apidb.load()
    typetable.register_classes(e["name"] for e in entries if e["kind"] in ("class", "struct")
                               and e["name"] not in overrides.IGNORED_CLASSES)
    typetable.register_result_structs(e["name"] for e in entries if e.get("fields"))
    records = [e for e in entries if e["kind"] in ("class", "struct")
               and e["name"] not in overrides.IGNORED_CLASSES]
    cpp, jl = emit(records)
    emitutil.write_cpp("types", cpp)
    emitutil.write_jl("types", jl)
    print(f"types: {len(records)} records -> out/cpp/types.cpp, out/julia/types.jl")


if __name__ == "__main__":
    main()
