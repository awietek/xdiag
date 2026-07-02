# SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0
#
# DMRG reference ground-state energy for the Boson "long" block test.
#
# Model: Bose-Hubbard chain, OPEN boundary conditions
#
#     H = -t * sum_{i} (b^dag_i b_{i+1} + b^dag_{i+1} b_i)
#         + (U/2) * sum_{i} n_i (n_i - 1)
#
# Local Hilbert space truncated at maxOcc bosons per site (dim d = maxOcc).
# Fixed total particle number: exactly `Nb` bosons.
#
# NOTE: with a fixed particle number a chemical potential -mu*sum n_i only
# shifts the energy by -mu*Nb, so it is dropped here. If you re-add it, remember
# to subtract mu*Nb before comparing.
#
# Matching xdiag block: Boson(N, d, Nb). The Boson block implements only "Id"
# and "Matrix" operators, so the C++ test builds the hopping / interaction from
# Op("Matrix", ...) local matrices in the d-dimensional on-site basis.
#
# Run:  julia --project dmrg_boson.jl

using ITensors
using ITensorMPS
using Printf

let
    N      = 65
    t      = 1.0
    U      = 1.0
    maxOcc = 4          # local boson dimension d
    Nb     = 2          # total number of bosons (low, fixed)

    sites = siteinds("Boson", N; dim=maxOcc, conserve_qns=true)

    ops = OpSum()
    for i in 1:N-1
        ops += -t, "Adag", i, "A", i + 1
        ops += -t, "Adag", i + 1, "A", i
    end
    for i in 1:N
        # U/2 n(n-1) = U/2 n^2 - U/2 n
        ops += U / 2.0, "N", i, "N", i
        ops += -U / 2.0, "N", i
    end
    H = MPO(ops, sites)

    # Initial product state with exactly Nb bosons, spread one-per-site.
    state = fill("0", N)
    for k in 1:Nb
        state[clamp(round(Int, (k - 0.5) * N / Nb), 1, N)] = "1"
    end
    psi0 = MPS(sites, state)

    nsweeps = 200
    maxdim  = [50, 100, 200, 400, 800, 1600]
    cutoff  = [1e-14]
    noise   = [1e-6, 1e-7, 1e-8, 0.0]

    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, noise)

    @printf("\nMODEL   Boson Bose-Hubbard OBC\n")
    @printf("N=%d  t=%.3f  U=%.3f  maxOcc=%d  Nb=%d\n", N, t, U, maxOcc, Nb)
    @printf("ENERGY  %.14f\n", energy)
end
