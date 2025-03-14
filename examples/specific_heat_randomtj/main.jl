using XDiag
using Random, Distributions
using LinearAlgebra
using Printf
using HDF5


function get_H!(H::OpSum, Nsites::Int, thopp::Float64, Jhopp::Float64)
    #Variance of the hopping and interaction parameter:
    thopp = 1.0
    Jhopp = 1.0
    Jdist = Normal(0.0, Jhopp^2) # Normal distribution
    tdist = Normal(0.0, thopp^2 / sqrt(2.0)) # thopp

    for i in 1:Nsites
        for j in (i+1):Nsites
            H["J_{$i}_{$j}"] = rand(Jdist) / sqrt(Nsites)
            H["t_{$i}_{$j}"] = (rand(tdist) + 1im * rand(tdist)) / sqrt(Nsites)
        end
    end
end


function main()
    Random.seed!(123) # seed this

    Nsamples = 200
    # define Hamiltonian
    Nsites = 8 # total number of sites
    nup = 2 # total number of spin-up fermions
    ndn = 2 # total number of spin-down fermions
    block = tJ(Nsites, nup, ndn) # tj model sites

    Temp = LinRange(0.01, 0.5, 64) # temperature array
    C = zeros(Float64, (length(Temp), Nsamples)) # array to store specific heat

    #build Hamiltonian
    H = OpSum()
    for i in 1:Nsites
        for j in (i+1):Nsites
            H += "J_{$i}_{$j}" * Op("SdotS", [i, j])
            H += "t_{$i}_{$j}" * Op("Hop", [i, j])
        end
    end

    for i in 1:Nsamples
        get_H!(H, Nsites, 1.0, 1.0) #generate random Hamiltonian
        Hmat = matrix(H, block)
        evals = eigvals(Hermitian(Hmat))

        for (j, t) in enumerate(Temp)
            exp_eval = exp.(-evals ./ t)
            Z = sum(exp_eval)  # Partition function
            energy = sum(evals .* exp_eval) / Z
            C[j, i] = (1.0 / t^2) * (sum(evals .* evals .* exp_eval) / Z - energy^2)
        end
    end

    filename = @sprintf("randomtj.Nsites.%d.outfile.h5", Nsites)

    h5open(filename, "w") do file
        write(file, "Temperature", collect(Temp))
        write(file, "Specific_Heat", C)
    end
end

main()