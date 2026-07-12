#!/bin/bash
project_root=/Users/awietek/Research/Software/xdiag
julia_prefix_path=`julia -e 'using CxxWrap; print(CxxWrap.prefix_path())'`
echo Julia prefix path: $julia_prefix_path

current_dir=`pwd`
cd $project_root

cmake -S . -B build -D XDIAG_JULIA_WRAPPER=On -D CMAKE_PREFIX_PATH=$julia_prefix_path -D CMAKE_CXX_COMPILER=clang++ -D CMAKE_EXPORT_COMPILE_COMMANDS=On

cd $current_dir
