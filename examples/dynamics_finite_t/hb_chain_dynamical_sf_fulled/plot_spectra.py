#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pydiag as yd
import pydiag.ensemble as yde

n_sites = 12
qs = [3] #range(n_sites//2+1)    # momenta q

temperatures = [0.001] # [0.001, 0.2, 0.5, 1.0]

n_omegas = 100           # number of frequencies 
max_omega = 4.0          # maximum frequency
eta = 0.05               # broadening factor

omegas = np.linspace(0, max_omega, n_omegas)
min_weight = 1e-10

# define ensemble of quantum numbers with degeneracies (qn, deg)
nups = [(nup, 2) if nup != n_sites // 2 else (nup, 1) for nup in range(n_sites//2+1)]
qs = [(q, 1) if q == 0 or q == n_sites // 2 else (q, 2) for q in range(n_sites//2+1)]
ensemble = yde.Ensemble(nups, qs)

directory = "outfiles/N.{}".format(n_sites)
regex = "outfile.N.{}.nup.(.*).k.(.*).seed.{}.iters.{}.h5".format(n_sites, seed, iters)

for q in qs:
    print("q =", q)
    S_of_q_tag = "S_of_q_{}_eig".format(q)
    data = yd.read_h5_data(directory, regex, tags=["eigenvalues", S_of_q_tag])
    eigs = yde.Array(ensemble, data, tag="eigenvalues")
    S_of_q_eig = yde.Array(ensemble, data, tag=S_of_q_tag)
    S_of_q_eig_dag = yde.conj(yde.transpose(A_tilde))

    for T in temperatures:
        print("T =", T)
        beta = 1 / T
        boltzmann = np.exp(-beta * eigs)
        partition = np.sum(boltzmann)

        poles = -yde.subtract_outer(eigs, eigs).flatten()
        weights = (2 * np.pi  / partition) * \
                  yde.einsum("a,ab,ba->ab", boltzmann,
                            S_of_q_eig_dag, S_of_q_eig).flatten()

        # Throw away poles with negligible weight
        poles = poles[np.abs(weights) > max_weight]
        weights = np.real(weights[np.abs(weights) > min_weight])

        # broaden and plot
        diffs = np.subtract.outer(omegas, poles)
        gaussians = np.exp(-(diffs / (2*eta))**2) / (eta * np.sqrt(2*np.pi))
        spectrum = gaussians @ weights
        plt.plot(omegas, spectrum, label=r"$q={}, T={:.3f}$".format(q, T))

plt.xlabel(r"$\omega$")
plt.ylabel(r"$S(q,\omega)$")
plt.title(r"$N={}$".format(n_sites))
plt.legend()
plt.show()
