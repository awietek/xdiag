#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

nsites = 16
nup = nsites // 2

eta = 0.05
n_omegas = 100
max_omega = 4.0
omegas = np.linspace(0, max_omega, n_omegas)

# Compute eigenvalues and eigenvectors of T-matrix
eigs = []
evecs = []
norms = []
e0s = []
for q in range(nsites):
    alphas = np.loadtxt("outfiles/alphas.N.{}.nup.{}.q.{}.txt".format(nsites, nup, q))
    betas = np.loadtxt("outfiles/betas.N.{}.nup.{}.q.{}.txt".format(nsites, nup, q))
    norm = np.loadtxt("outfiles/norm.N.{}.nup.{}.q.{}.txt".format(nsites, nup, q))

    # compute poles and weights from tridiagonal matrix
    tmat = np.diag(alphas[:50]) + np.diag(betas[:49], k=1) + np.diag(betas[:49], k=-1)
    eig, evec = np.linalg.eigh(tmat)
    eigs.append(eig)
    evecs.append(evec)
    norms.append(norm)
    e0s.append(eig[0])
e0 = min(e0s)
    
# compute broadened spectrum with Gaussian broadening eta
spectra = np.zeros((nsites, n_omegas))
for q in range(nsites):
    poles = eigs[q] - e0
    weights = (norms[q]**2) * (evecs[q][0, :]**2)
    diffs = np.subtract.outer(omegas, poles)
    gaussians = np.exp(-(diffs / (2*eta))**2) / (eta * np.sqrt(2*np.pi))
    spectra[q, :] = gaussians @ weights

p = plt.imshow(spectra.T, origin='lower', extent=[0, nsites, 0, max_omega],
           interpolation=None, aspect='auto')
plt.xlabel(r"$q$")
plt.ylabel(r"$\omega$")
plt.colorbar(p)
plt.show()
