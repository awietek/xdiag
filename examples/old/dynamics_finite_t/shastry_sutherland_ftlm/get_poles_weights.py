#!/usr/bin/env python
import os
import numpy as np
import h5py
import pydiag as yd
import pydiag.ensemble as yde
from collections import OrderedDict

nsites = 16
q = "M.C1.A"

J=0.63
Jd=1.00

temperature = 0.500

n_omegas = 100           # number of frequencies 
max_omega = 4.0          # maximum frequency
eta = 0.05               # broadening factor

niter = 200


seeds = range(1, 51, 1)
cutoffs = [1e-8, 1e-10, 1e-12]
niter_evals = [100, 150, 200]  

# seeds = [1]
# cutoffs = [1e-4]
# niter_evals = [200]  



# define ensemble of quantum numbers with degeneracies (qn, deg)
nups = [(nup, 1) if nup == nsites // 2 else (nup, 2) for nup in range(nsites//2+1)]
ks = ["Gamma.C1.A", "M.C1.A", "X0.C1.A", "X1.C1.A"]

# nups = [(1, 1)]
# ks = ["X1.C1.A"]
ensemble = yde.Ensemble(nups, ks)

beta = 1.0 / temperature

print("T = {:.4f}".format(temperature))
print("q = ", q)

for seed in seeds:
    for cutoff in cutoffs:
        for niter_eval in niter_evals:
            print("@ seed: {}, cutoff: {}, niter: {}".format(seed, cutoff, niter_eval))
            directory = "/data/condmat/awietek/Research/Software/xdiag/examples/dynamics_finite_t/shastry_sutherland_ftlm/outfiles/shastry.{}.HB.J.Jd.fsl/J.{:.2f}.Jd.{:.2f}/q.{}/seed.{}/".format(nsites, J, Jd, q, seed)
            regex = "outfile.shastry.{}.HB.J.Jd.fsl.J.{:.2f}.Jd.{:.2f}.q.{}.seed.{}.nup.(.*).k.(.*).niter.{}.h5".format(nsites, J, Jd, q, seed, niter)

            print("Reading data")
            data = yd.read_h5_data(directory, regex, tags=["AlphasT", "BetasT", "DimK", "AlphasS", "BetasS", 'A'])
            dims = yde.Scalar(ensemble, data, tag="DimK")
            T_tridiag = yde.TriDiag(ensemble, data, diag_tag="AlphasT", offdiag_tag="BetasT")
            T = T_tridiag.asarray()
            S_tridiag = yde.TriDiag(ensemble, data, diag_tag="AlphasS", offdiag_tag="BetasS")
            S = S_tridiag.asarray()
            A = yde.Array(ensemble, data, tag="A")

            # Adjust matrices
            print("Adjusting matrix size")
            for block, deg in ensemble:
                # resize
                T.array[block] = T.array[block][:niter_eval, :niter_eval]
                S.array[block] = S.array[block][:niter_eval, :niter_eval]
                A.array[block] = A.array[block][:niter_eval, :niter_eval]
                T_tridiag.diag[block] = T_tridiag.diag[block][:niter_eval]
                T_tridiag.offdiag[block] = T_tridiag.offdiag[block][:niter_eval-1]
                S_tridiag.diag[block] = S_tridiag.diag[block][:niter_eval]
                S_tridiag.offdiag[block] = S_tridiag.offdiag[block][:niter_eval-1]

            # Compute totel e0
            e0s = T_tridiag.eig0()
            e0 = e0s.min()
            print("e0:", e0)

            print("computing d, Q, e, R")
            d, Q = T_tridiag.eig()
            e, R = S_tridiag.eig()
            Q_dag = yde.transpose(Q)
            R_dag = yde.transpose(R)
            
            # Adjust matrices
            print("Computing T - e0, S - e0")
            for block, deg in ensemble:
                T_size = T.array[block].shape[0]
                T.array[block] -= e0 * np.eye(T_size)
                S_size = S.array[block].shape[0]
                S.array[block] -= e0 * np.eye(S_size)
                d.array[block] -= e0
                e.array[block] -= e0            

                
            print("computing boltzmann factors b, b_tilde, partition")
            eT = yde.expm(-beta/2 * T)
            b = eT[:,0]
            b_tilde = yde.dot(Q_dag, b)
            b_tilde_dag = yde.transpose(b_tilde)

            # compute partition
            partitions = yde.dot(b, b) * dims
            partition = 0.
            for block, deg in ensemble:
                partition += partitions.array[block]
            
            print("computing A_tilde")
            A_tilde_tmp = yde.dot(A, Q)
            A_tilde = yde.dot(R_dag, A_tilde_tmp)
            A_tilde_dag = yde.conj(yde.transpose(A_tilde))
    
            print("computing weights")
            weights = 2 * np.pi * dims * yde.einsum("a,ab,bc,c->abc", \
                                                    b_tilde_dag, A_tilde_dag,\
                                                    A_tilde, b_tilde).flatten() 
    
            print("computing poles")
            poles = yde.add_outer(-0.5 * d, yde.add_outer(e, -0.5 * d)).flatten()

            # security check
            for block, deg in ensemble:
                if weights.array[block].shape != poles.array[block].shape:
                    raise ValueError("size of weights and poles is not the same")

            print("Concatenating poles and weights")
            all_weights = weights.concatenate() / partition
            all_poles = poles.concatenate()    
            print("npoles:         ", len(all_poles))

            if cutoff != None:
                print("Truncating poles and weights")
                inds = np.abs(all_weights) > cutoff
                all_poles = all_poles[inds]
                all_weights = all_weights[inds]
                print("npoles: ({:.1e})".format(cutoff), len(all_poles))    

            print("Sorting poles and weights")
            inds = all_poles.argsort()
            all_weights_sorted = all_weights[inds]
            all_poles_sorted = all_poles[inds]

            print("imaginary part norm:", np.linalg.norm(np.imag(all_weights_sorted)))
            # block= ('1', 'X1.C1.A')
            # print(d.array[block])
            # print(e.array[block])
            # for p, w in zip(all_poles_sorted, all_weights_sorted):
            #     print("p: {:.16f} w: {}".format(p, w))
            # print()

            # print(b_tilde_dag.array[block])
            # print(b_tilde.array[block])
            
            # for a_idx, dd1 in enumerate(d.array[block]):
            #     for b_idx, ee in enumerate(e.array[block]):
            #         for c_idx, dd2 in enumerate(d.array[block]):
            #             p = ee - (dd1 + dd2) / 2
            #             f1 = b_tilde_dag.array[block][a_idx]
            #             f2 = A_tilde_dag.array[block][a_idx, b_idx]
            #             f3 = A_tilde.array[block][b_idx, c_idx]
            #             f4 = b_tilde.array[block][c_idx]
            #             w = f1 * f2 * f3 * f4
            #             print("a: {} b: {} c: {}".format(a_idx, b_idx, c_idx))
            #             print("p: {:.6f} w: {:.3e}".format(p, w))
            #             print("f1: {:.3e} f2: {:.3e} f3: {:.3e} f4: {:.3e}".format(f1, f2, f3, f4))
            #             print()
                        
            directory = "/data/condmat/awietek/Research/Software/xdiag/examples/dynamics_finite_t/shastry_sutherland_ftlm/poles_weights/shastry.{}.HB.J.Jd.fsl/J.{:.2f}.Jd.{:.2f}/q.{}/T.{:.4f}/niter.{}.cutoff.{}".format(nsites, J, Jd, q, temperature, niter_eval, cutoff)
            filename = "poles.weights.shastry.{}.HB.J.Jd.fsl.J.{:.2f}.Jd.{:.2f}.q.{}.T.{:.4f}.niter.{}.cutoff.{}.seed.{}.h5".format(nsites, J, Jd, q, temperature, niter_eval, cutoff, seed)
            if not os.path.exists(directory):
                os.makedirs(directory)
    
            filename = os.path.join(directory, filename)
            with h5py.File(filename, "w") as fl:
                fl["poles"] = all_poles_sorted

                # if complex, weights always occur in conjugate pairs
                # -> total imaginary part is zero anyway, can be ignored
                fl["weights"] = np.real(all_weights_sorted)  
                fl["partition"] = partition
