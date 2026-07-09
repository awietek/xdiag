#!/bin/bash

project_root="/Users/awietek/Research/Software/xdiag"
generated_dir="/Users/awietek/Research/Software/xdiag/misc/wrapgen2/generated/julia"

mkdir -p $generated_dir

./code/emit_types_julia.py > $generated_dir/types.jl

modules=`ls -d $project_root/xdiag/*/ | sed "s|$project_root/xdiag/||" | sed "s|/||"`

for mod in ${modules[@]}; do
    echo $mod
    ./code/emit_submodule_julia.py $mod > $generated_dir/$mod.jl
done

# install armadillo wrapper
cp $project_root/misc/wrapgen2/static/armadillo.jl $generated_dir
./code/emit_xdiag_jl.py > $generated_dir/XDiag.jl

# Format
julia -e "using JuliaFormatter; format(\"$generated_dir\")"
