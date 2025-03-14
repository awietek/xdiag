using LinearAlgebra, LinearAlgebra.LAPACK
using XDiag
using HDF5
using Printf
const Nsites = 20 # Total number of sites in the Shastry-Sutherland model:
const Rtqp = 2 # Number of random vectors to be used in the algorithm

function Tqp_mag_sector(Obser::Matrix{Float64}, Temp::LinRange{Float64,Int64}, H::OpSum)

    block = Spinhalf(Nsites) #Create spin-1/2 block with conservation of Sz

    for k in 1:Rtqp #Perform the calculation for each random vector
        res = eigvals_lanczos(H, block, neigvals=1, precision=1e-12, max_iterations=150, deflation_tol=1e-7, random_seed=k) # Perform the Lanczos interation starting from a random vector;

        d = length(res.alphas)
        vecs = Matrix{Float64}(I, d, d)
        #This part can be done in the post-processing
        (eigs, vecs) = LAPACK.stev!('V', res.alphas, res.betas)

        for k in 1:length(Temp)
            psi1 = exp.(-(eigs .- eigs[1]) ./ (2.0 .* Temp[k])) .* conj(vecs[1, :]) # compute psi_1

            Obser[k, 1] += dot(psi1, psi1)# get norm factor for the partition function

            Obser[k, 2] += dot(psi1, eigs .* psi1) # measure energy

            Obser[k, 3] += dot(psi1, eigs .* (eigs .* psi1))# measure energy square
        end
    end
end

function main()


    #Reads the file with the interactions
    fl = FileToml("shastry_sutherland_L_5_W_4.toml")
    ops = read_opsum(fl, "Interactions")
    #Defines exchange coupling strenght -- This ratio of couplings corresponds to the dimmer phase.
    ops["Jd"] = 1.0
    ops["J"] = 0.630

    # Linear Array with target temperatures
    Temp = LinRange(0.01, 0.35, 64)

    # Array to store the average value of the observable 
    Obser = zeros(Float64, (length(Temp), 3))

    Tqp_mag_sector(Obser, Temp, ops)

    filename = @sprintf("shastry_sutherland_L_5_W_4.h5")
    h5open(filename, "w") do file
        write(file, "Temp", collect(Temp))
        write(file, "Observable", Obser)
    end
end

main()