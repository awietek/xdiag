# SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0
#
# DMRG reference ground-state energy for the Electron "long" block test.
#
# Model: single-band Hubbard chain, OPEN boundary conditions
#
#     H = -t * sum_{i,sigma} (c^dag_{i sigma} c_{i+1 sigma} + h.c.)
#         + U * sum_{i} n_{i up} n_{i dn}
#
# Fixed particle numbers: exactly `nup` up and `ndn` down electrons.
#
# Matching xdiag OpSum (built in the C++ test):
#     for i in 0..N-2:
#         ops += t * Op("Hop", {i, i+1})  # xdiag Hop = -sum_sigma(c^dag c+h.c.)
#     ops += U * Op("HubbardU")           # HubbardU = sum_i n_{i up} n_{i dn}
#     block = Electron(N, nup, ndn)
#
# Run:  julia --project dmrg_electron.jl

using ITensors
using ITensorMPS
using Printf

let
    N   = 65
    t   = 1.0
    U   = 4.0
    nup = 1
    ndn = 1

    sites = siteinds("Electron", N; conserve_qns=true)

    ops = OpSum()
    for i in 1:N-1
        ops += -t, "Cdagup", i, "Cup", i + 1
        ops += -t, "Cdagup", i + 1, "Cup", i
        ops += -t, "Cdagdn", i, "Cdn", i + 1
        ops += -t, "Cdagdn", i + 1, "Cdn", i
    end
    for i in 1:N
        ops += U, "Nupdn", i
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
    noise   = [1e-6, 1e-7, 1e-8, 0.0]

    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, noise)

    @printf("\nMODEL   Electron Hubbard OBC\n")
    @printf("N=%d  t=%.3f  U=%.3f  nup=%d  ndn=%d\n", N, t, U, nup, ndn)
    @printf("ENERGY  %.14f\n", energy)
end
