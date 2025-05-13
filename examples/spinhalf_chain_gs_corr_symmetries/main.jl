using XDiag

function main()

    # periodic N-site Heisenberg antiferromagnet
    N = 14
    H = OpSum()
    for i in 1:N
        H += "J" * Op("SdotS", [i, mod1(i+1, N)])
    end
    H["J"] = 1.0 

    # define spin correlator at "half-chain" distance
    corr_op = OpSum()
    corr_op += Op("SdotS", [1, N ÷ 2])

    # define shift by one site
    T = Permutation(circshift(1:N, -1))
    
    # define cyclic group
    C_N = PermutationGroup([T^k for k in 0:(N-1)])
    
    # define irreps from character tables, labelled by momentum (2pi/N ×) k
    character_table = [ [C_N_character(N, k, p) for p in 0:(N-1)] for k in 0:(N-1)]
    irreps = [ Representation(C_N, character_table[k]) for k in 1:N ]   

    # check all translation symmetry blocks
    e0 = Inf
    psi0 = nothing
    for k in eachindex(irreps)
        block = Spinhalf(N, irreps[k]) 
        e0k, psik = eig0(H, block)
        if  e0k < e0
            e0 = e0k
            psi0 = psik
        end
    end

    # now the operator must be symmetrized!
    corr = inner(symmetrize(corr_op, C_N), psi0)

    println("Ground state energy: $e0")
    println("Ground state correlator: $corr")
    println("State vector length = ", size(psi0))

end

function C_N_character(N::Int, k::Int, p::Int)
    exp(im * 2 * pi * k * p / N)
end

main()




