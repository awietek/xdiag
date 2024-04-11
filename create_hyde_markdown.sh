#!/bin/bash

xdiagdir="/Users/awietek/Research/Software/xdiag"
includes="-I$xdiagdir -I/opt/homebrew/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX13.3.sdk/usr/include"

sources=(`find $xdiagdir/xdiag -name *.hpp -not -path "*old*"`)
blacklist=combinatorics,detail,utils,terms,mpi,basis,random,omp,mpi,electron,spinhalf,tj,bitops,symmetries,operators,numeric,variant

target_dir="autoreference"
mkdir -p $target_dir

for src in ${sources[@]}; do
    echo Compiling $src
    hyde --hyde-update -hyde-yaml-dir=$target_dir --namespace-blacklist=$blacklist $src -- -std=c++17 -c $includes
done
