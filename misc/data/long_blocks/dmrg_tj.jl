# SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0
#
# DMRG reference ground-state energy for the tJ "long" block test.
#
# Model: t-J chain, OPEN boundary conditions (no double occupancy)
#
#     H = -t * sum_{i,sigma} (c^dag_{i sigma} c_{i+1 sigma} + h.c.)
#         + J * sum_{i} ( S_i . S_{i+1} - n_i n_{i+1} / 4 )
#
# Fixed particle numbers: exactly `nup` up and `ndn` down electrons.
#
# Matching xdiag OpSum (built in the C++ test):
#     for i in 0..N-2:
#         ops += t * Op("Hop",     {i, i+1})  # xdiag Hop = -sum_sigma(c^dag c+h.c.)
#         ops += J * Op("tJSdotS", {i, i+1})  # tJSdotS = S.S - n_i n_j / 4
#     block = tJ(N, nup, ndn)
#
# Run:  julia --project dmrg_tj.jl

using ITensors
using ITensorMPS
using Printf

let
    N   = 65
    t   = 1.0
    J   = 1.0
    nup = 1
    ndn = 1

    sites = siteinds("tJ", N; conserve_qns=true)

    ops = OpSum()
    for i in 1:N-1
        # hopping, both spins, both directions
        ops += -t, "Cdagup", i, "Cup", i + 1
        ops += -t, "Cdagup", i + 1, "Cup", i
        ops += -t, "Cdagdn", i, "Cdn", i + 1
        ops += -t, "Cdagdn", i + 1, "Cdn", i
        # J (S_i.S_{i+1} - n_i n_{i+1} / 4)
        ops += J, "Sz", i, "Sz", i + 1
        ops += J / 2.0, "S+", i, "S-", i + 1
        ops += J / 2.0, "S-", i, "S+", i + 1
        ops += -J / 4.0, "Ntot", i, "Ntot", i + 1
    end
    H = MPO(ops, sites)

    # Initial product state: nup "Up" and ndn "Dn" on distinct, spread sites.
    state = fill("Emp", N)
    positions = [clamp(round(Int, (k - 0.5) * N / (nup + ndn)), 1, N)
                 for k in 1:(nup + ndn)]
    for k in 1:nup
        state[positions[k]] = "Up"
    end
    for k in 1:ndn
        state[positions[nup + k]] = "Dn"
    end
    psi0 = MPS(sites, state)

    nsweeps = 200
    maxdim  = [50, 100, 200, 400, 800, 1600]
    cutoff  = [1e-14]
    noise   = [1e-6, 1e-6,1e-6,1e-6,1e-6, 1e-6,1e-6, 1e-6,1e-6,1e-6,
               1e-7, 1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,
               1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8]

    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, noise)

    @printf("\nMODEL   tJ nearest-neighbour OBC\n")
    @printf("N=%d  t=%.3f  J=%.3f  nup=%d  ndn=%d\n", N, t, J, nup, ndn)
    @printf("ENERGY  %.14f\n", energy)
end
