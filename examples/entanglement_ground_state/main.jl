using LinearAlgebra
using XDiag
using HDF5
using Printf

function push_ProductState!(ps1, psA, psB)
    for sub in psA
        push!(ps1, String(sub))
    end
    for sub in psB
        push!(ps1, String(sub))
    end
end

function Create_RDM(RDM, vstate, Block, BlockA, BlockB) #Function to construct the RDM
    for (i, p1) in enumerate(BlockA)
        for (j, p2) in enumerate(BlockA)
            for (k, p3) in enumerate(BlockB)
                p1aux = ProductState()
                p2aux = ProductState()
                push_ProductState!(p1aux, p1, p3)
                push_ProductState!(p2aux, p2, p3)
                id1 = index(Block, p1aux)
                id2 = index(Block, p2aux)
                RDM[j, i] += vstate[id2] * conj(vstate[id1])
            end
        end
    end
end

function getVnE(RDM) #Function to compute the entanglement entropy.
    vne = 0
	p = eigvals(Hermitian(RDM))#FULL ED in the RDM
    for i in p
        if i>1e-14
            vne -= i*log(i)
        end
    end
	return vne
end

function main()
    N = 16 # chain length
    ops = OpSum()
    block = Spinhalf(N) # Spin1/2 Block full
    Delta = 1.0
    for i in 1:N
        ops += Op("Exchange", [i, mod1(i + 1, N)])
        ops += Delta*Op("SzSz", [i, mod1(i + 1, N)])
    end

	e0, gs = eig0(ops, block)# get GS with Lanczos	
    vstate = vector(gs)
	vne = zeros(Float64,(length(1:div(N,2)),2))
    for Na in 1:div(N,2)
		Nb = N - Na
		blockA = Spinhalf(Na) # Spin1/2 Block full
		blockB = Spinhalf(Nb) # Spin1/2 Block full
		RDM = zeros(ComplexF64, (dim(blockA), dim(blockA)))
		Create_RDM(RDM, vstate, block, blockA, blockB)
		vne[Na,1] = Na
		vne[Na,2] = getVnE(RDM)
	end
	filename = @sprintf("EE_XXX_model.Nsites.%d.Delta.%d.outfile.h5", N,Delta)
    h5open(filename, "w") do file
        write(file, "vne", vne)
    end	
end

main()