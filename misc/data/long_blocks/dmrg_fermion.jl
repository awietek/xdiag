# SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0
#
# DMRG reference ground-state energy for the Fermion "long" block test.
#
# Model: spinless fermions, nearest-neighbour hopping + interaction, OPEN chain
#
#     H = -t * sum_{i} (c^dag_i c_{i+1} + c^dag_{i+1} c_i)
#         + V * sum_{i} n_i n_{i+1}
#
# Fixed particle number: exactly `Nf` fermions.
#
# Matching xdiag OpSum (built in the C++ test):
#     for i in 0..N-2:
#         ops += t * Op("Hop", {i, i+1})   # xdiag Hop = -(c^dag c + h.c.), so
#                                          # coupling +t gives physical -t hopping
#         ops += V * Op("NN",  {i, i+1})   # NN = n_i n_j
#     block = Fermion(N, Nf)
#
# Run:  julia --project dmrg_fermion.jl

using ITensors
using ITensorMPS
using Printf

let
    N  = 65
    t  = 1.0
    V  = 1.0
    Nf = 2          # number of fermions (low, fixed)

    sites = siteinds("Fermion", N; conserve_qns=true)

    ops = OpSum()
    for i in 1:N-1
        ops += -t, "Cdag", i, "C", i + 1
        ops += -t, "Cdag", i + 1, "C", i
        ops += V, "N", i, "N", i + 1
    end
    H = MPO(ops, sites)

    # Initial product state with exactly Nf occupied sites, spread evenly.
    state = fill("Emp", N)
    for k in 1:Nf
        state[clamp(round(Int, (k - 0.5) * N / Nf), 1, N)] = "Occ"
    end
    psi0 = MPS(sites, state)

    nsweeps = 200
    maxdim  = [50, 100, 200, 400, 800, 1600]
    cutoff  = [1e-14]
    noise   = [1e-6, 1e-7, 1e-8, 0.0]

    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, noise)

    @printf("\nMODEL   Fermion spinless t-V OBC\n")
    @printf("N=%d  t=%.3f  V=%.3f  Nf=%d\n", N, t, V, Nf)
    @printf("ENERGY  %.14f\n", energy)
end
