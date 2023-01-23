#!/bin/bash

n_sites=12

n_ups=`seq 0 1 $(( n_sites / 2 ))`
ks=`seq 0 1 $(( n_sites / 2))` 

# n_ups=(2) #`seq 0 1 $(( n_sites / 2 ))`
# ks=(0)    #`seq 0 1 $(( n_sites / 2))` 


for n_up in ${n_ups[@]}; do
for k in ${ks[@]}; do

    ./main $n_sites $n_up $k

done
done
