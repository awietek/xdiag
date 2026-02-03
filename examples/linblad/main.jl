using XDiag
using LinearAlgebra
using Printf
using HDF5

Sx = [0.0 0.5; 0.5 0.0]
Sy = [0.0 -0.5*1im; 0.5*1im 0.0]


function main()
    L = 6 # total number of sites

    gamma = 2.0 # dissipation strength

    Nsites_liouville_space = 2 * L ## the effective Hilbert Space is doubled!
    
    nup = Nsites_liouville_space ÷ 2
    
    block = Spinhalf(Nsites_liouville_space, nup)

    #build Hamiltonian
    H = OpSum()
    for i in 1:L
        H += -1.0*1im*Op("Matrix", [i, mod1(i+1, L)],kron(Sx,Sx)) ## normal spin chain
        H += -1.0*1im*Op("Matrix", [i, mod1(i+1, L)],kron(Sy,Sy))

        H += 1.0im*Op("Matrix", [L+i, L+mod1(i+1, L)],kron(Sx,Sx)) ## normal spin chain
        H += 1.0im*Op("Matrix", [L+i, L+mod1(i+1, L)],kron(Sy,Sy))

        H += gamma*Op("SzSz", [i, L+i]) ## recycling term
    end
        
    
    Hmat = matrix(H, block)

    evals = eigvals(Hmat)    

    evals .-= gamma*L/4

    filename = @sprintf("./heisenberg_bulkd_dephasing.Nsites.%d.gamma.%.3f.nup_sector.%d.outfile.h5", L,gamma,nup)

    h5open(filename, "w") do file
        write(file, "real_part_eigenvalues", real.(evals))
        write(file, "imag_part_eigenvalues", imag.(evals))
    end
end

main()