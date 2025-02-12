#!/usr/bin/env python
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
import pydiag as yd
import pydiag.ensemble as yde
from collections import OrderedDict

nsites = 20
J=1.00
Jd=1.00

seeds = range(1, 21, 1)
temperatures = np.linspace(0.01, 1.0, 100)

# define ensemble of quantum numbers with degeneracies (qn, deg)
nups = [(nup, 1) if nup == nsites // 2 else (nup, 2) for nup in range(nsites//2+1)]
if nsites == 16:
    ks = [("Gamma.D4.A1", 1), ("Gamma.D4.A2", 1), ("Gamma.D4.B1", 1),
          ("Gamma.D4.B2", 1), ("Gamma.D4.E", 2), ("M.D4.A1", 1),
          ("M.D4.A2", 1), ("M.D4.B1", 1), ("M.D4.B2", 1), ("M.D4.E", 2),
          ("Sigma.D1.A", 4), ("Sigma.D1.B", 4), ("X.D2.A1", 2),
          ("X.D2.A2", 2), ("X.D2.B1", 2), ("X.D2.B2", 2)]

if nsites == 20:
    ks = [("Gamma.C4.A", 1), ("Gamma.C4.B", 1), ("Gamma.C4.Ea", 1),
          ("Gamma.C4.Eb", 1), ("M.C4.A", 1), ("M.C4.B", 1), ("M.C4.Ea", 1),
          ("M.C4.Eb", 1), ("None0.C1.A", 4), ("None1.C1.A", 4)]

ensemble = yde.Ensemble(nups, ks)

for seed in seeds:
    print("seed:", seed)
    directory = "outfiles/seed.{}".format(seed)
    regex = "outfile.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.nup.(.*).k.(.*).seed.{}.h5".format(nsites, J, Jd, seed)
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
      
        eT = yde.expm(-beta * tmat_e0 / 2)  ## Division by 2 bnecause of TPQ state
        b = eT[:,0]

        partition = yde.sum(dims * yde.dot(b, b), degeneracies=True)
        energy = yde.sum(dims * yde.dot(b, yde.dot(tmat_arr, b)), degeneracies=True)
        energy2 = yde.sum(dims * yde.dot(yde.dot(tmat_arr, b), yde.dot(tmat_arr, b)), degeneracies=True)
        print(partition, energy / partition, energy2 / partition)
        output[Tidx, 0] = T
        output[Tidx, 1] = partition
        output[Tidx, 2] = energy / partition
        output[Tidx, 3] = energy2 / partition

        
    np.savetxt("moments/moments.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.seed.{}.txt".format(nsites, J, Jd, seed), output, 
               header = "{:<24} {:<24} {:<24} {:<24}".format("T", "Z", "E", "E2"))
