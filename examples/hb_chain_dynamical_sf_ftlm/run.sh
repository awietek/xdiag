#!/bin/bash

n_sites=16
n_ups=`seq 0 1 $n_sites`
ks=`seq 0 1 $(( n_sites - 1))` 
seeds=`seq 1 1 10`

n_ups=(0)
ks=(1) 
seeds=(1)

iters=200

for seed in ${seeds[@]}; do
for n_up in ${n_ups[@]}; do
for k in ${ks[@]}; do

    ./main $n_sites $n_up $k $seed $iters    

done
done
done
