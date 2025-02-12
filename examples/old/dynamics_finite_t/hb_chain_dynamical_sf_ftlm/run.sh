#!/bin/bash

nsites=16

nups=`seq 0 1 $nsites`
ks=`seq 0 1 $(( nsites - 1))` 
seeds=`seq 1 1 10`

# nups=(0)
# ks=(1) 
# seeds=(1)

iters=200

for seed in ${seeds[@]}; do
for nup in ${nups[@]}; do
for k in ${ks[@]}; do

    ./main $nsites $nup $k $seed $iters    

done
done
done
