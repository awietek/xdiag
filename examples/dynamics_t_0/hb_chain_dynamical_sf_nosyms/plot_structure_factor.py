#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

nsites = 12
nup = nsites // 2

eta = 0.05
n_omegas = 100
max_omega = 4.0
omegas = np.linspace(0, max_omega, n_omegas)

# compute broadened spectrum
spectra = np.zeros((nsites, n_omegas))
for q in range(nsites):
    alphas = np.loadtxt("outfiles/alphas.N.{}.nup.{}.q.{}.txt".format(nsites, nup, q))
    betas = np.loadtxt("outfiles/betas.N.{}.nup.{}.q.{}.txt".format(nsites, nup, q))
    norm = np.loadtxt("outfiles/norm.N.{}.nup.{}.q.{}.txt".format(nsites, nup, q))

    # compute poles and weights from tridiagonal matrix
    tmat = np.diag(alphas[:50]) + np.diag(betas[:49], k=1) + np.diag(betas[:49], k=-1)
    eigs, evecs = np.linalg.eigh(tmat)
    e0 = eigs[0]
    poles = eigs - e0
    weights = (norm**2) * (evecs[0, :]**2)
    
    # compute broadened spectrum with Gaussian broadening eta
    diffs = np.subtract.outer(omegas, poles)
    gaussians = np.exp(-(diffs / (2*eta))**2) / (eta * np.sqrt(2*np.pi))
    spectra[q, :] = gaussians @ weights

p = plt.imshow(spectra.T, origin='lower', extent=[0, nsites, 0, max_omega],
           interpolation=None, aspect='auto')
plt.xlabel(r"$q$")
plt.ylabel(r"$\omega$")
plt.colorbar(p)
plt.show()
