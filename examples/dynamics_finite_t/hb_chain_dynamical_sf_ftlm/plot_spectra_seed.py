#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

n_sites = 20
q = 3
T = 0.1
seeds = [1, 2, 3, 4]
etas = [0.05]
prec = "1e-12"

n_omegas = 100
max_omega = 4.0
omegas = np.linspace(0, max_omega, n_omegas)

for eta in etas:
    for seed in seeds:
        filename = "data/spectrum.N.{}.q.{}.T.{:.3f}.seed.{}.eta.{}.prec.{}.npy".format(n_sites, q, T, seed, eta, prec)
        spectrum = np.load(filename)
        plt.plot(omegas, np.imag(spectrum), label=r"$\eta={}$, seed {}".format(eta, seed))
plt.legend()
plt.title(r"$N={}, T={:.3f}, q={}$".format(n_sites, T, q))
plt.show()
        
