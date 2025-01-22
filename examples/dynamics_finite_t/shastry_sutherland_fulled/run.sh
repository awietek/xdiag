#!/bin/bash

nsites=16
nups=`seq 0 1 $(( nsites / 2 ))`
ks=("Gamma.C1.A" "M.C1.A" "X0.C1.A" "X1.C1.A") 
qs=("M.C1.A") 

J=0.63
Jd=1.00


for nup in ${nups[@]}; do
for k in ${ks[@]}; do
for q in ${qs[@]}; do

    ./build/main $nsites $nup $k $q $J $Jd

done
done
done
