#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pydiag as yd
import pydiag.ensemble as yde
from collections import OrderedDict

nsites = 12
qs = [3] #range(nsites//2+1)    # momenta q

temperatures = [0.001, 0.2, 0.5, 1.0]

n_omegas = 100           # number of frequencies 
max_omega = 4.0          # maximum frequency
eta = 0.05               # broadening factor

omegas = np.linspace(-max_omega, max_omega, n_omegas)
min_weight = 1e-8

# define ensemble of quantum numbers with degeneracies (qn, deg)
nups = [(nup, 1) if nup == nsites // 2 else (nup, 2) for nup in range(nsites//2+1)]
ks = [(k, 1) if k == 0 or k == nsites // 2 else (k, 2) for k in range(nsites//2+1)]
ensemble = yde.Ensemble(nups, ks)

directory = "outfiles/N.{}".format(nsites)
regex = "outfile.N.{}.nup.(.*).k.(.*).h5".format(nsites)

for q in qs:
    print("q =", q)
    S_of_q_tag = "S_of_q_{}_eig".format(q)
    data = yd.read_h5_data(directory, regex, tags=["eigenvalues", S_of_q_tag])
    eigs = yde.Array(ensemble, data, tag="eigenvalues").flatten()
    e0 = eigs.min().min()
    eigs -= e0

    # define eigenvalues with momentum p + q
    data_q = OrderedDict()
    for (nup, p), deg in ensemble:
        pq = int(p) + int(q)
        if (pq > nsites // 2):
            pq = nsites - (int(p) + int(q))
        pq = str(pq)
        # print(p, q, pq)
        # print(eigs.array.keys())
        data_q[(nup, p)] = eigs.array[(nup, pq)]
    eigs_q = yde.Array(ensemble, data_q)
    
    S_of_q_eig = yde.Array(ensemble, data, tag=S_of_q_tag)
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
