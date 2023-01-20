#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

n_sites = 12
qs = [3] #range(n_sites//2+1)    # momenta q

n_up = n_sites // 2

eta = 0.05
n_omegas = 100
max_omega = 4.0
omegas = np.linspace(0, max_omega, n_omegas)

# compute broadened spectrum
spectra = np.zeros((n_sites, n_omegas))
for q in qs:
    alphas = np.loadtxt("outfiles/alphas.N.{}.nup.{}.q.{}.txt".format(n_sites, n_up, q))
    betas = np.loadtxt("outfiles/betas.N.{}.nup.{}.q.{}.txt".format(n_sites, n_up, q))
    norm = np.loadtxt("outfiles/norm.N.{}.nup.{}.q.{}.txt".format(n_sites, n_up, q))

    # compute poles and weights from tridiagonal matrix
    tmat = np.diag(alphas) + np.diag(betas[:-1], k=1) + np.diag(betas[:-1], k=-1)
    eigs, evecs = np.linalg.eigh(tmat)
    e0 = eigs[0]
    poles = eigs - e0
    weights = 2* np.pi * (norm**2) * (evecs[0, :]**2)

    # compute broadened spectrum with Gaussian broadening eta
    diffs = np.subtract.outer(omegas, poles)
    gaussians = np.exp(-(diffs / (2*eta))**2) / (eta * np.sqrt(2*np.pi))
    spectrum = gaussians @ weights
    plt.plot(omegas, spectrum, label=r"$q={}$".format(q))

plt.xlabel(r"$\omega$")
plt.ylabel(r"$S(q,\omega)$")
plt.title(r"$N={}$".format(n_sites))
plt.legend()
plt.show()
