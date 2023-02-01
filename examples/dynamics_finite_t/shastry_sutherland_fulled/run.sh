#!/bin/bash

n_sites=16
n_ups=`seq 0 1 $(( n_sites / 2 ))`
ks=("Gamma.C1.A" "M.C1.A" "X0.C1.A" "X1.C1.A") 
qs=("M.C1.A") 

J=0.63
Jd=1.00


for n_up in ${n_ups[@]}; do
for k in ${ks[@]}; do
for q in ${qs[@]}; do

    ./main $n_sites $n_up $k $q $J $Jd

done
done
done
