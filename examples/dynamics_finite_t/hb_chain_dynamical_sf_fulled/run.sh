#!/bin/bash

nsites=12

nups=`seq 0 1 $(( nsites / 2 ))`
ks=`seq 0 1 $(( nsites / 2))` 

# nups=(2) #`seq 0 1 $(( nsites / 2 ))`
# ks=(0)    #`seq 0 1 $(( nsites / 2))` 


for nup in ${nups[@]}; do
for k in ${ks[@]}; do

    ./build/main $nsites $nup $k

done
done
