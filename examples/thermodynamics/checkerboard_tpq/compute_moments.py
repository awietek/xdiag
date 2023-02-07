#!/usr/bin/env python
import os
import numpy as np
import matplotlib.pyplot as plt
import pydiag as yd
import pydiag.ensemble as yde
from collections import OrderedDict

n_sites = 20
J=1.00
Jd=1.00

seeds = range(1, 21, 1)
temperatures = np.linspace(0.01, 1.0, 2)

seeds = [1]
temperatures = [1.0]



# define ensemble of quantum numbers with degeneracies (qn, deg)
nups = [(nup, 1) if nup == n_sites // 2 else (nup, 2) for nup in range(n_sites//2+1)]
if n_sites == 16:
    ks = [("Gamma.D4.A1", 1), ("Gamma.D4.A2", 1), ("Gamma.D4.B1", 1),
          ("Gamma.D4.B2", 1), ("Gamma.D4.E", 2), ("M.D4.A1", 1),
          ("M.D4.A2", 1), ("M.D4.B1", 1), ("M.D4.B2", 1), ("M.D4.E", 2),
          ("Sigma.D1.A", 4), ("Sigma.D1.B", 4), ("X.D2.A1", 2),
          ("X.D2.A2", 2), ("X.D2.B1", 2), ("X.D2.B2", 2)]
    
if n_sites == 20:
    ks = [("Gamma.C4.A", 1), ("Gamma.C4.B", 1), ("Gamma.C4.Ea", 1),
          ("Gamma.C4.Eb", 1), ("M.C4.A", 1), ("M.C4.B", 1), ("M.C4.Ea", 1),
          ("M.C4.Eb", 1), ("None0.C1.A", 4), ("None1.C1.A", 4)]
    
ensemble = yde.Ensemble(nups, ks)

for seed in seeds:
    print("seed:", seed)
    directory = "outfiles/seed.{}".format(seed)
    regex = "outfile.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.nup.(.*).k.(.*).seed.{}.h5".format(n_sites, J, Jd, seed)
    data = yd.read_h5_data(directory, regex)
    tmat = yde.TriDiag(ensemble, data, diag_tag="Alphas", offdiag_tag="Betas")
    eigs = yde.Array(ensemble, data, tag="Eigenvalues").flatten()
    dims = yde.Scalar(ensemble, data, tag="Dimension")
    e0 = eigs.min().min()

    tmat_arr = tmat.asarray()
    output = np.zeros((len(temperatures), 4))
    
    for Tidx, T in enumerate(temperatures):
        print("T =", T)
        beta = 1 / T
        tmat_e0 = tmat_arr - e0 * yde.eye_like(tmat_arr)

        for block, d in ensemble:
            print(block)
            # print(tmat_arr.array[block])
            # print(tmat_e0.array[block])
            print(eigs.array[block])
            print(np.linalg.eigvalsh(tmat_arr.array[block]))
            print(np.linalg.eigvalsh(tmat_e0.array[block]))
            print()
        
        eT = yde.expm(-beta * tmat_e0)
        b = eT[:,0]
        
        partition = yde.sum(dims * yde.dot(b, b), degeneracies=True)
        energy = yde.sum(dims * yde.dot(b, yde.dot(tmat_arr, b)), degeneracies=True)
        energy2 = yde.sum(dims * yde.dot(yde.dot(tmat_arr, b), yde.dot(tmat_arr, b)), degeneracies=True)
        print(partition, energy / partition, energy2 / partition)
        output[Tidx, 0] = T
        output[Tidx, 1] = partition
        output[Tidx, 2] = energy / partition
        output[Tidx, 3] = energy2 / partition

    np.savetxt("moments/moments.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.seed.{}.txt".format(n_sites, J, Jd, seed), output, 
               header = "{:<24} {:<24} {:<24} {:<24}".format("T", "Z", "E", "E2"))
