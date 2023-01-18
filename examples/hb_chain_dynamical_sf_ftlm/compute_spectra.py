#!/usr/bin/env python
import numpy as np
import pickle
import pydiag as yd

n_sites = 16
iters = 200

seed = 1
q=1

beta = 1.0

directory = "outfiles/N.{}".format(n_sites)
regex = "outfile.N.{}.nup.(.*).k.(.*).seed.{}.iters.{}.h5".format(n_sites, seed, iters)

data = yd.read_h5_data(directory, regex, tags=["T_alphas", "T_betas"])
nups = [(nup, 2) if nup != n_sites // 2 else (nup, 1) for nup in range(3)]
qs = [(q, 1) if q == 0 or q == n_sites // 2 else (q, 2) for q in range(n_sites//2+1)]


# nups = [(nup, 2) if nup != n_sites // 2 else (nup, 1) for nup in range(n_sites//2+1)]
# qs = [(q, 1) if q == 0 or q == n_sites // 2 else (q, 2) for q in range(n_sites//2+1)]

print("computing e0")
ensemble = yd.Ensemble(nups, qs)
T_tridiag = yd.EnsembleTriDiag(ensemble, data, diag_tag="T_alphas", offdiag_tag="T_betas")
e0s = T_tridiag.eig0()
e0 = e0s.min()
print(e0)

print("computing b")
T = T_tridiag.asarray()
eT = yd.ensemble_expm( -beta * (T - e0))
b = eT[:,0]


print("computing b_tilde")
d, Q = T_tridiag.eig()
Q_dag = yd.ensemble_transpose(Q)
b_tilde = yd.ensemble_dot(Q_dag, b)
b_tilde_dag = yd.ensemble_transpose(b)

print("computing R")
S_diag_tag = "q_{}/S_alphas".format(q)
S_offdiag_tag = "q_{}/S_betas".format(q)
A_tag = "q_{}/A".format(q)
data = yd.read_h5_data(directory, regex,
                       tags=[S_diag_tag, S_offdiag_tag, A_tag])

S_tridiag = yd.EnsembleTriDiag(ensemble, data, diag_tag=S_diag_tag,
                               offdiag_tag=S_offdiag_tag,)
S = S_tridiag.asarray()
e, R = S_tridiag.eig()

print("computing A_tilde")
A = yd.EnsembleMatrix(ensemble, data, tag=A_tag)
A_tilde = yd.ensemble_dot(Q_dag, yd.ensemble_dot(A, R))
A_tilde_dag = yd.ensemble_conj(yd.ensemble_transpose(A_tilde))

    
# print("min", ens_e0.min())
# for block, deg in ensemble:
#     print("blocks", block, deg)


    
# pkl = "data.pkl"
# with open(pkl, "wb") as f: # "wb" because we want to write in binary mode
#     pickle.dump(data, f)


# with open(pkl, "rb") as f:
#     data2 = pickle.load(f)
#     # print(data2)
