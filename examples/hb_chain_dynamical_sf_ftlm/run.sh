#!/bin/bash

n_sites=20
n_ups=`seq 0 1 $n_sites`
seeds=`seq 1 1 5`

iters=200

for n_up in ${n_ups[@]}; do
    for seed in ${seeds[@]}; do
    ./main $n_sites $n_up $seed $iters    
    done
done


