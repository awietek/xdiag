#!/usr/bin/env python
import numpy as np
import h5py
import pydiag as yd
import pydiag.ensemble as yde

nsites = 16
iters = 200

seed = 1

beta = 1.0

directory = "outfiles/N.{}".format(nsites)
regex = "outfile.N.{}.nup.(.*).k.(.*).seed.{}.iters.{}.h5".format(nsites, seed, iters)

# nups = [(nup, 2) if nup != nsites // 2 else (nup, 1) for nup in range(3)]
# qs = [(q, 1) if q == 0 or q == nsites // 2 else (q, 2) for q in range(nsites//2+1)]
nups = [(nup, 2) if nup != nsites // 2 else (nup, 1) for nup in range(nsites//2+1)]
qs = [(q, 1) if q == 0 or q == nsites // 2 else (q, 2) for q in range(nsites//2+1)]
ensemble = yde.Ensemble(nups, qs)

print("Reading T")
data = yd.read_h5_data(directory, regex, tags=["T_alphas", "T_betas", "dim"])

print("computing e0")
T_tridiag = yde.TriDiag(ensemble, data, diag_tag="T_alphas", offdiag_tag="T_betas")
e0s = T_tridiag.eig0()

dims = yde.Scalar(ensemble, data, tag="dim")
for block, e0, dim in zip(e0s.keys(), e0s.values(), dims.values()):
    print(block, e0, dim)

e0 = e0s.min()
print(e0)

print("computing b")
T = T_tridiag.asarray()
eT = yde.expm( -beta * (T - e0))
b = eT[:,0]


print("computing b_tilde")
d, Q = T_tridiag.eig()
Q_dag = yde.transpose(Q)
b_tilde = yde.dot(Q_dag, b)
b_tilde_dag = yde.transpose(b)

q=1

print("Reading S and A")
S_diag_tag = "q_{}/S_alphas".format(q)
S_offdiag_tag = "q_{}/S_betas".format(q)
A_tag = "q_{}/A".format(q)
data = yd.read_h5_data(directory, regex,
                       tags=[S_diag_tag, S_offdiag_tag, A_tag])
S_tridiag = yde.TriDiag(ensemble, data, diag_tag=S_diag_tag,
                        offdiag_tag=S_offdiag_tag,)
S = S_tridiag.asarray()

print("computing e and R")
e, R = S_tridiag.eig()

print("computing A_tilde")
A = yde.Array(ensemble, data, tag=A_tag)
A_tilde = yde.dot(Q_dag, yde.dot(A, R))
A_tilde_dag = yde.conj(yde.transpose(A_tilde))

print("computing weights")
weights = 2 * np.pi * dims * yde.einsum("a,ab,bc,c->abc", \
                                        b_tilde_dag, A_tilde,\
                                        A_tilde_dag, b_tilde).flatten() 
print(weights[])

print("computing poles")
poles = yde.add_outer(-0.5 * e, yde.add_outer(d, -0.5 * e)).flatten()
print(poles)
# a = np.ones((2,))
# b = np.ones((2,2))
# c = np.ones((2,2))
# d = np.ones((2,))
# e = np.einsum("a,ab,bc,c->abc", a,b,c,d)
# print(e)

# print("min", ens_e0.min())
# for block, deg in ensemble:
#     print("blocks", block, deg)


    
# pkl = "data.pkl"
# with open(pkl, "wb") as f: # "wb" because we want to write in binary mode
#     pickle.dump(data, f)


# with open(pkl, "rb") as f:
#     data2 = pickle.load(f)
#     # print(data2)
