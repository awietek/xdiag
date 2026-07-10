#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
# Regenerate the whole wrapper: model -> types -> submodules -> top-level wiring.
set -e
cd "$(dirname "$0")"

python3 model.py                       # xdiag headers -> api/*.toml

mkdir -p out/cpp out/julia
cp static/*.hpp static/*.cpp out/cpp/                   # static C++ bridges
cp static/*.jl out/julia/                               # static Julia bridges

python3 emit_types.py                  # api/*.toml -> out/{cpp,julia}/types.*
python3 emit_submodule.py              # api/*.toml -> out/{cpp,julia}/<module>.* + coverage.txt
python3 emit_module.py                 # wire it all -> xdiagjl.{hpp,cpp}, XDiag.jl

# Create cmake source file
cmake_source_file=out/cpp/sources.cmake
echo "Writing CMake source file -> $cmake_source_file"
echo "set(XDIAG_JULIA_SOURCES" > $cmake_source_file
for fl in out/cpp/*.cpp; do
    echo $fl | sed 's|out/cpp|    src|' >> $cmake_source_file
done
echo ")" >> $cmake_source_file

clang-format -i out/cpp/*.cpp out/cpp/*.hpp
julia -e 'using JuliaFormatter; format("out/julia")'
echo "done. review out/coverage.txt for what was skipped."
