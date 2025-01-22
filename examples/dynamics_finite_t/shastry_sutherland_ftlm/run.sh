#!/bin/bash

nsites=16
nups=`seq 0 1 $(( nsites / 2 ))`
ks=("Gamma.C1.A" "M.C1.A" "X0.C1.A" "X1.C1.A") 
qs=("M.C1.A") 

latfile="shastry.$nsites.HB.J.Jd.fsl.toml"
J=0.63
Jd=1.00

niter=500
seeds=`seq 1 1 50`

for q in ${qs[@]}; do
for seed in ${seeds[@]}; do
for nup in ${nups[@]}; do
for k in ${ks[@]}; do

    outdir="/data/condmat/awietek/Research/Software/xdiag/examples/dynamics_finite_t/shastry_sutherland_ftlm/outfiles/shastry.$nsites.HB.J.Jd.fsl/J.$J.Jd.$Jd/q.$q/seed.$seed"
    mkdir -p $outdir
    outfile="$outdir/outfile.shastry.$nsites.HB.J.Jd.fsl.J.$J.Jd.$Jd.q.$q.seed.$seed.nup.$nup.k.$k.niter.$niter.h5"

    dumpdir="/scratch/awietek/Research/Software/xdiag/examples/dynamics_finite_t/shastry_sutherland_ftlm/dumpfiles/shastry.$nsites.HB.J.Jd.fsl/J.$J.Jd.$Jd/q.$q/seed.$seed/nup.$nup.k.$k.niter.$niter"
    mkdir -p $dumpdir    
    
    ./build/main \
	--latfile $latfile \
	--J $J \
	--Jd $Jd \
	--nup $nup \
	--k $k \
	--q $q \
	--seed $seed \
	--niter $niter \
	--outfile $outfile \
	--dumpdir $dumpdir 

    rm -rf $dumpdir    
done
done
done
done
