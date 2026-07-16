# SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0
#
# DMRG reference ground-state energy for the Spinhalf "long" block test.
#
# Model: spin-1/2 Heisenberg chain, OPEN boundary conditions
#
#     H = J * sum_{i=0}^{N-2} S_i . S_{i+1}
#       = J * sum_i [ Sz_i Sz_{i+1} + 1/2 (S+_i S-_{i+1} + S-_i S+_{i+1}) ]
#
# Fixed magnetization sector: exactly `nup` up-spins (total Sz = nup - N/2).
#
# Matching xdiag OpSum (built in the C++ test):
#     for i in 0..N-2:  ops += J * Op("SdotS", {i, i+1})
#     block = Spinhalf(N, nup)
#
# Run:  julia --project dmrg_spinhalf.jl
# Record the printed "ENERGY" value into the C++ test.

using ITensors
using ITensorMPS
using Printf

let
    N   = 65
    J   = 1.0
    nup = 2          # number of up-spins (low, fixed -> small sector)

    sites = siteinds("S=1/2", N; conserve_qns=true)

    ops = OpSum()
    for i in 1:N-1
        ops += J,       "Sz", i, "Sz", i + 1
        ops += J / 2.0, "S+", i, "S-", i + 1
        ops += J / 2.0, "S-", i, "S+", i + 1
    end
    H = MPO(ops, sites)

    # Initial product state with exactly nup up-spins, spread evenly to help
    # convergence into the target magnetization sector.
    state = fill("Dn", N)
    for k in 1:nup
        state[clamp(round(Int, (k - 0.5) * N / nup), 1, N)] = "Up"
    end
    psi0 = MPS(sites, state)

    nsweeps = 200
    maxdim  = [50, 100, 200, 400, 800, 1600]
    cutoff  = [1e-14]
    noise   = [1e-6, 1e-7, 1e-8, 0.0]

    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, noise)

    @printf("\nMODEL   Spinhalf Heisenberg OBC\n")
    @printf("N=%d  J=%.3f  nup=%d\n", N, J, nup)
    @printf("ENERGY  %.14f\n", energy)
end
