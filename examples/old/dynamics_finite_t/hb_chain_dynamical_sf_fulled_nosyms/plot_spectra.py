#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py

nsites = 12
qs = [3] #range(nsites//2+1)    # momenta q

temperatures = [0.001, 0.2, 0.5, 1.0]

n_omegas = 100           # number of frequencies 
max_omega = 4.0          # maximum frequency
eta = 0.05               # broadening factor

omegas = np.linspace(-max_omega, max_omega, n_omegas)
min_weight = 1e-8

outfile = "outfiles/N.{}/outfile.N.{}.h5".format(nsites, nsites)

with h5py.File(outfile, 'r') as fl:
    eigs = fl["eigenvalues"][0]
eigs = eigs - eigs[0]

for q in qs:
    print("q =", q)
    with h5py.File(outfile, 'r') as fl:
        S_of_q_tmp = fl["S_of_q_{}_eig".format(q)][:]
        S_of_q_eig = S_of_q_tmp['real']  + 1j*S_of_q_tmp['imag']
    S_of_q_eig_dag = np.conj(np.transpose(S_of_q_eig))

    for T in temperatures:
        print("T =", T)
        beta = 1 / T
        boltzmann = np.exp(-beta * eigs)
        partition = np.sum(boltzmann)
        print(partition)
        poles = -np.subtract.outer(eigs, eigs).flatten()
        weights = (2 * np.pi  / partition) * \
                  np.einsum("a,ab,ba->ab", boltzmann,
                            S_of_q_eig_dag, S_of_q_eig).flatten()

        # Throw away poles with negligible weight
        poles = poles[np.abs(weights) > min_weight]
        weights = np.real(weights[np.abs(weights) > min_weight])

        print(np.sort(poles))
        
        # broaden and plot
        diffs = np.subtract.outer(omegas, poles)
        gaussians = np.exp(-(diffs / (2*eta))**2) / (eta * np.sqrt(2*np.pi))
        spectrum = gaussians @ weights
        plt.plot(omegas, spectrum, label=r"$q={}, T={:.3f}$".format(q, T))
        print(spectrum)
        
plt.xlabel(r"$\omega$")
plt.ylabel(r"$S(q,\omega)$")
plt.title(r"$N={}$".format(nsites))
plt.legend()
plt.show()
