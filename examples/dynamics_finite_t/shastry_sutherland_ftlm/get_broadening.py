#!/usr/bin/env python
import os
import numpy as np
import h5py
import pydiag as yd
import pydiag.ensemble as yde
from collections import OrderedDict
import matplotlib.pyplot as plt
from scipy.stats import sem

nsites = 16
q = "M.C1.A"

J=0.63
Jd=1.00

temperature = 0.5

n_omegas = 100           # number of frequencies 
max_omega = 4.0          # maximum frequency
eta = 0.01               # broadening factor

niter = 200

# seeds = range(1, 51, 1)
# cutoffs = [1e-8, 1e-10, 1e-12]
# niters = [100, 150, 200]  

seeds = range(1, 21, 1)
cutoffs = [1e-8]
niters = [200]  

# define ensemble of quantum numbers with degeneracies (qn, deg)
nups = [(nup, 1) if nup == nsites // 2 else (nup, 2) for nup in range(nsites//2+1)]
ks = ["Gamma.C1.A", "M.C1.A", "X0.C1.A", "X1.C1.A"]


ensemble = yde.Ensemble(nups, ks)

beta = 1.0 / temperature

print("T = {:.4f}".format(temperature))
print("q = ", q)

n_omegas = 100           # number of frequencies 
max_omega = 4.0          # maximum frequency
eta = 0.1               # broadening factor
omegas = np.linspace(0, max_omega, n_omegas)


for cutoff in cutoffs:
    for niter in niters:
        spectra = np.zeros((len(seeds), len(omegas)))

        for seed_idx, seed in enumerate(seeds):

            directory = "/data/condmat/awietek/Research/Software/xdiag/examples/dynamics_finite_t/shastry_sutherland_ftlm/poles_weights/shastry.{}.HB.J.Jd.fsl/J.{:.2f}.Jd.{:.2f}/q.{}/T.{:.4f}/niter.{}.cutoff.{}".format(nsites, J, Jd, q, temperature, niter, cutoff)
            filename = "poles.weights.shastry.{}.HB.J.Jd.fsl.J.{:.2f}.Jd.{:.2f}.q.{}.T.{:.4f}.niter.{}.cutoff.{}.seed.{}.h5".format(nsites, J, Jd, q, temperature, niter, cutoff, seed)
            if not os.path.exists(directory):
                os.makedirs(directory)
    
            filename = os.path.join(directory, filename)
            with h5py.File(filename, "r") as fl:
                poles = fl["poles"][:]
                weights = fl["weights"][:]    
                partition = fl["partition"]

            diffs = np.subtract.outer(omegas, poles)
            gaussians = np.exp(-(diffs / (2*eta))**2) / (eta * np.sqrt(2*np.pi))
            spectrum = gaussians @ weights * nsites**2
            spectra[seed_idx, :] = spectrum
            plt.plot(omegas, spectrum, label=r"$seed={}, cutoff={}, niter={}$".format(seed, cutoff, niter))

        spectrum_mean = np.mean(spectra, axis=0)
        spectrum_err = sem(spectra, axis=0)
        plt.errorbar(omegas, spectrum_mean, spectrum_err, c="k", lw=3, capsize=2)

plt.legend()
plt.show()
