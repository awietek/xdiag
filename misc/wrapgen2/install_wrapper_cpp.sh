#!/bin/bash
project_root="/Users/awietek/Research/Software/xdiag"
generated_dir="/Users/awietek/Research/Software/xdiag/misc/wrapgen2/generated/cpp"
artifact_hash="44722d868723dcb2366a9e90f02b34885cd80aba"

rm -rf $project_root/julia/src
mkdir -p $project_root/julia/src
cp $generated_dir/* $project_root/julia/src

origin=`pwd`
cd $project_root

# cmake -S . -B build -D XDIAG_JULIA_WRAPPER=On -D CMAKE_PREFIX_PATH=/Users/awietek/.julia/artifacts/44722d868723dcb2366a9e90f02b34885cd80aba -D CMAKE_CXX_COMPILER=clang++

cmake -S . -B build -D XDIAG_JULIA_WRAPPER=On -D CMAKE_PREFIX_PATH=/Users/awietek/.julia/artifacts/$arifact_hash -D CMAKE_CXX_COMPILER=clang++ -D CMAKE_EXPORT_COMPILE_COMMANDS=On
cmake --build build -j10
cmake --install build

# removes Julia's potentially faulty precompilation cache
# see troubleshooting.md for details
rm -rf ~/.julia/compiled/*/XDiag    

cd $origin
