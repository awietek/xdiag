#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pydiag as yd
import pydiag.ensemble as yde
from collections import OrderedDict

n_sitess = [16, 20]
J=1.00
Jd=1.00

temperatures = np.linspace(0.0001, 10.0, 100)

# define ensemble of quantum numbers with degeneracies (qn, deg)

for n_sites in n_sitess:

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

    directory = "outfiles/".format(n_sites)
    regex = "outfile.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.nup.(.*).k.(.*).h5".format(n_sites, J, Jd)

    data = yd.read_h5_data(directory, regex, tags=["Eigenvalues"])
    eigs = yde.Array(ensemble, data, tag="Eigenvalues").flatten()
    e0 = eigs.min().min()
    eigs -= e0

    specheats = []
    for T in temperatures:
        print("T =", T)
        beta = 1 / T
        boltzmann = yde.exp(-beta * eigs)
        partition = yde.sum(boltzmann, degeneracies=True)
        energy = yde.sum(eigs * boltzmann, degeneracies=True)
        energy2 = yde.sum(eigs * eigs* boltzmann, degeneracies=True)
        specheat = (energy2/partition - (energy/partition)**2) / T
        specheats.append(specheat)

    plt.plot(temperatures, specheats, label=r"$N={}$".format(n_sites))

plt.legend()
plt.xlabel(r"$T$")
plt.ylabel(r"$C(T)$")
plt.title(r"$J={:.2f}, J_d={:.2f}$".format(J, Jd))
plt.legend()
plt.show()
