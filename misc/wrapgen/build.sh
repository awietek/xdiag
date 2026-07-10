#!/bin/bash

project_root=/Users/awietek/Research/Software/xdiag
origin_source_dir=$project_root/misc/wrapgen/out/cpp
target_source_dir=$project_root/julia/src
prefix_path=$(julia -e 'using CxxWrap; print(CxxWrap.prefix_path())')
echo Julia prefix path: $prefix_path

rm $target_source_dir/*
rsync -a --checksum $origin_source_dir/ $target_source_dir/

current_dir=`pwd`
cd $project_root

cmake -S . -B build -D XDIAG_JULIA_WRAPPER=On -D CMAKE_PREFIX_PATH=$prefix_path -D CMAKE_CXX_COMPILER=clang++ -D CMAKE_EXPORT_COMPILE_COMMANDS=On
#-D CMAKE_BUILD_TYPE=Debug

cmake --build build -j10

cd $current_dir
