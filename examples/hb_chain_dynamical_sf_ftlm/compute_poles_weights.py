#!/usr/bin/env python
import numpy as np
import h5py

n_sites = 20
seed = 1

T = 0.1
eigs = dict()
for nup in range(n_sites + 1):
    outfile = "outfiles/N.{}.nup.{}/seed.{}/outfile.h5".format(n_sites, nup, seed)
    with h5py.File(outfile, 'r') as fl:
        alphas = fl["alphas"][0]
        betas = fl["betas"][0]
        tmat = np.diag(alphas) + np.diag(betas[:-1], k=1) + np.diag(betas[:-1], k=-1)
        eigs[nup] = np.linalg.eigvalsh(tmat)
    
    e0 = min([eig[0] for eig in eigs.values()])
    print("nup: {} e0: {}", nup, e0)
    

    
