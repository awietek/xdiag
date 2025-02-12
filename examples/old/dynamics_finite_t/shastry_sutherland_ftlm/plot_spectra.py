#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pydiag as yd
import pydiag.ensemble as yde
from collections import OrderedDict

nsites = 16
q = "M.C1.A"

J=0.63
Jd=1.00

T = 0.1

n_omegas = 100           # number of frequencies 
max_omega = 4.0          # maximum frequency
eta = 0.05               # broadening factor

omegas = np.linspace(0, max_omega, n_omegas)
min_weight = 1e-8

# define ensemble of quantum numbers with degeneracies (qn, deg)
nups = [(nup, 1) if nup == nsites // 2 else (nup, 2) for nup in range(nsites//2+1)]
ks = ["Gamma.C1.A", "M.C1.A", "X0.C1.A", "X1.C1.A"]
ensemble = yde.Ensemble(nups, ks)

directory = "outfiles/".format(nsites)
regex = "outfile.shastry.{}.J.{:.2f}.Jd.{:.2f}.nup.(.*).k.(.*).q.{}.h5".format(nsites, J, Jd, q)

print("q =", q)
data = yd.read_h5_data(directory, regex, tags=["EigenvaluesK", "EigenvaluesKQ", "SofQ"])
eigs = yde.Array(ensemble, data, tag="EigenvaluesK").flatten()
e0 = eigs.min().min()

print(eigs.min())
eigs -= e0
eigs_q = yde.Array(ensemble, data, tag="EigenvaluesKQ").flatten()
print(eigs_q.min())

eigs_q -= e0


S_of_q_eig = yde.Array(ensemble, data, tag="SofQ")
S_of_q_eig_dag = yde.conj(yde.transpose(S_of_q_eig))

for T in temperatures:
    print("T =", T)
    beta = 1 / T
    boltzmann = yde.exp(-beta * eigs)
    partition = yde.sum(boltzmann, degeneracies=True)
    print(partition)
    poles = -yde.subtract_outer(eigs, eigs_q).flatten()
    weights = (2 * np.pi  / partition) * \
              yde.einsum("a,ab,ba->ab", boltzmann,
                         S_of_q_eig_dag, S_of_q_eig).flatten()

    # Turn ensemble array into np.array
    poles = poles.concatenate(degeneracies=True)
    weights = weights.concatenate(degeneracies=True)

    # Throw away poles with negligible weight
    poles = poles[np.abs(weights) > min_weight]
    weights = np.real(weights[np.abs(weights) > min_weight])

    print(np.sort(poles))

    # broaden and plot
    diffs = np.subtract.outer(omegas, poles)
    gaussians = np.exp(-(diffs / (2*eta))**2) / (eta * np.sqrt(2*np.pi))
    spectrum = gaussians @ weights * nsites**2
    plt.plot(omegas, spectrum, label=r"$q={}, T={:.3f}$".format(q, T))
    print(spectrum)
        
plt.xlabel(r"$\omega$")
plt.ylabel(r"$S(q,\omega)$")
plt.title(r"$N={}$".format(nsites))
plt.legend()
plt.show()
