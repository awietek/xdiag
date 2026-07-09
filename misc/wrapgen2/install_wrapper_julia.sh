#!/bin/bash
project_root="/Users/awietek/Research/Software/xdiag"
generated_dir="/Users/awietek/Research/Software/xdiag/misc/wrapgen2/generated/julia"
target_dir="/Users/awietek/.julia/dev/XDiag/src"

cp $generated_dir/* $target_dir
cp static/armadillo.jl $target_dir
