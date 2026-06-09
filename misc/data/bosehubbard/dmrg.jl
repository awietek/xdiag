using ITensors
using ITensorMPS

let 
    N=4
    t=1
    U=8
    mu=0.5

    maxOcc = 5

    # set up Hamiltonian for a PBC chain
    ops = OpSum()
    for s in 1:N
        s1 = s
        s2 = mod1(s+1, N)
        ops += -t, "Adag", s1, "A", s2
        ops += -t, "Adag", s2, "A", s1
        ops += U/2.0, "N", s, "N", s
        ops += -U/2.0, "N", s
        # ops += -mu,"N",s
    end

    # This says that the sites of the lattice have local bosons,
    # with maximal dimension maxOcc
    sites = siteinds("Boson", N; dim=maxOcc)

    # This creates a matrix-product operator for the Hamiltonian
    H = MPO(ops, sites)

    # to fix particle number specify initial state
    # https://docs.itensor.org/ITensorMPS/stable/tutorials/QN_DMRG.html
    psi0 = random_mps(sites)


    # DMRG specific parameters
    nsweeps = 1000      # number of DMRG sweeps
    maxdim = [1000]    # maximal bond dimension
    cutoff = [1E-16]   # how much to truncate teh density matrix
    energy, psi = dmrg(H,psi0; nsweeps, maxdim, cutoff)    # this runs DMRG
    @show energy

    M = correlation_matrix(psi, "Adag", "A")
    display(M)
end
