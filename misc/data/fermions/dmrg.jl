using ITensors
using ITensorMPS

let 
    N=6
    t=1
    V=3.7
    mu=-0.2

    # set up Hamiltonian for a PBC chain
    ops = OpSum()
    for s in 1:N-1
        s1 = s
        s2 = mod1(s+1, N)
        ops += -t, "Cdag", s1, "C", s2
        ops += -t, "Cdag", s2, "C", s1
        ops += V, "N", s1, "N", s2
    end

    for s in 1:N
        ops += mu,"N",s
    end
    
    # This says that the sites of the lattice have local bosons,
    # with maximal dimension maxOcc
    sites = siteinds("Fermion", N)

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

    M = correlation_matrix(psi, "Cdag", "C")
    display(M)
end
