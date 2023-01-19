#!/usr/bin/env python
import numpy as np
from scipy.linalg import expm
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import h5py

n_sites = 20
n_up = n_sites // 2

n_omegas = 100
max_omega = 4.0
omegas = np.linspace(0, max_omega, n_omegas)
T = 0.1

seeds=range(1, 6, 1)
qs = range(n_sites)
beta = 1 / T

prec = 1e-12
etas = [0.01, 0.02, 0.05, 0.1] 

# read_data
for seed in seeds:
    print("seed", seed)

    alphas = dict()
    betas = dict()
    eigs = dict()
    evecs = dict()
    tmats = dict()
    S_of_q = dict()

    # Read data
    for nup in range(n_sites):
        outfile = "outfiles/N.{}.nup.{}/seed.{}/outfile.h5".format(n_sites, nup, seed)
        # print("reading", outfile)
        with h5py.File(outfile, 'r') as fl:
            alphas[nup] = fl["alphas"][0]
            betas[nup] = fl["betas"][0]
            tmats[nup] = np.diag(alphas[nup]) + np.diag(betas[nup][:-1], k=1) + \
                                                        np.diag(betas[nup][:-1], k=-1)
            eigs[nup], evecs[nup] = np.linalg.eigh(tmats[nup])
            for q in qs:
                S_of_q[(q, nup)] = np.matrix(fl["s_of_q_{}".format(q)][:]['real'] + 1j*fl["s_of_q_{}".format(q)][:]['imag'])

    e0 = min([e[0] for e in eigs.values()])
    print("e0", e0)

    for q in qs:    
        print("computing poles and weigts, q = ", q)
        poles = []
        weights = []
        for nup in range(n_sites):
            iters = len(alphas[nup])
            tmat = tmats[nup] - e0 * np.eye(iters)
            Q = evecs[nup]

            b = expm( - beta * tmats[nup] / 2)[:,0]
\            b_tilde = Q.T @ b
            A = S_of_q[(q, nup)]
            A_tilde = (Q.T @ A) @ Q
            B = S_of_q[(q, nup)].H
            B_tilde = Q.T @ B @ Q

            eig = eigs[nup] - e0
            pls = np.add.outer(-eig/2, np.add.outer(eig, -eig/2)).flatten()
            wgts = 2 * np.pi * np.einsum("a,ab,bc,c->abc", b_tilde, A_tilde, B_tilde, b).flatten()

            idces = np.abs(wgts)>prec
            print("nup: ", nup, ", n_weights:", len(wgts[idces]))
            
            poles += list(pls[idces])
            weights += list(wgts[idces])
        poles = np.array(poles)
        weights = np.array(weights)

        # Broaden using scipy's kernel density estimation
        diffs = np.subtract.outer(omegas, poles)
        for eta in etas:
            print("broadening with eta = ", eta)
            gaussians = np.exp(-(diffs / (2*eta))**2) / (eta * np.sqrt(2*np.pi))
            spectrum = gaussians @ weights
            filename = "data/spectrum.N.{}.q.{}.T.{:.3f}.seed.{}.eta.{}.prec.{}".format(n_sites, q, T, seed, eta, prec)
            np.save(filename, spectrum)


# p = plt.imshow(spectra.T, origin='lower', extent=[0, n_sites, 0, max_omega],
#            interpolation=None, aspect='auto')
# plt.xlabel(r"$q$")
# plt.ylabel(r"$\omega$")
# plt.colorbar(p)
# plt.show()
