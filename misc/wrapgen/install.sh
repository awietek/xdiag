#!/bin/bash
project_root="/Users/awietek/Research/Software/xdiag"

origin_julia_dir=$project_root/misc/wrapgen/out/julia
target_julia_dir="/Users/awietek/.julia/dev/XDiag/src"

# Install cpp library
current_dir=`pwd`
cd $project_root
cmake --install build
cd $current_dir

# Install julia library
rm -r $target_julia_dir/*
cp $origin_julia_dir/* $target_julia_dir

# removes Julia's potentially faulty precompilation cache
# see troubleshooting.md for details
rm -rf ~/.julia/compiled/*/XDiag    
