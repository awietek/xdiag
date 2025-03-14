using XDiag
using LinearAlgebra
using Printf
using HDF5

function main()
    Nsites = 24 ## Number of Sites 
    alpha = 1.0 ## Exponent power-law decay
    ops = OpSum() ## create OpSum

    for i in 1:Nsites
        for j in (i+1):Nsites
        J = -1.0 / (sqrt(Nsites) * (abs(i - j))^(alpha))
        ops += J * Op("Exchange", [i, j])
        end
    end

    energies = Vector{Float64}[] ## Collect the energies
    
    Diagonalize in each magnetization sector using FullED
    for nup in 0:Nsites
        block = Spinhalf(Nsites, nup)
        H = matrix(ops, block)
        eig = eigvals(Hermitian(H))
        for e0 in eig
            push!(energies, [nup, e0])
        end
    end

    #Diagonalize in each magnetization sector using Lanczos to get the first few eigenvalues
    # for nup in 0:Nsites
    #     block = Spinhalf(Nsites, nup)
    #     r = eigvals_lanczos(ops, block, neigvals = 1);
    #     for e0 in r.eigenvalues
    #         push!(energies, [nup, e0])
    #     end
    # end
    
    filename = @sprintf("energies_tos_XXmodel.Nsites.%d.alpha.%d.outfile.h5", Nsites, alpha)
    h5open(filename, "w") do file
        write(file, "energies", hcat(energies...))
    end
end

main()