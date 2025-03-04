using XDiag

function main()

    # periodic N-site Heisenberg antiferromagnet
    N = 14
    H = OpSum()
    for i in 1:N
        H += "J" * Op("SdotS", [i, mod1(i+1, N)])
    end
    H["J"] = 1.0 

    # compute ground state correlator of this operator
    corr_op = OpSum()
    corr_op += Op("SdotS", [1, N ÷ 2])

    # run different implementations
    gs_correlator_simple(N, H, corr_op)

    gs_correlator_Sz_sym(N, H, corr_op)
    
    gs_correlator_translation_sym(N, H, corr_op)

end

# treat Hilbert space as a whole (memory consuming)
function gs_correlator_simple(N::Int, H::OpSum, corr_op::OpSum)
    full_block = Spinhalf(N) 
    lanczos_res = eigs_lanczos(H, full_block)
    e0 = lanczos_res.eigenvalues[1]
    psi0 = lanczos_res.eigenvectors    
    corr = inner(corr_op, psi0)
    
    println("Ground state energy: $e0")
    println("Ground state correlator: $corr")
    println("State vector length = ", size(psi0))
end

# use Sz_tot conservation (less memory consuming) 
function gs_correlator_Sz_sym(N::Int, H::OpSum, corr_op::OpSum)
    e0 = Inf
    psi0 = nothing
    Nup_vals = 0:N
    # check blocks with fixed number of up-spins
    for nup in Nup_vals
        sym_block = Spinhalf(N, nup) 
        lanczos_result = eigs_lanczos(H, sym_block)
        gs = lanczos_result.eigenvalues[1]
        if  gs < e0
            e0 = gs
            psi0 = lanczos_result.eigenvectors 
        end
    end    
    corr = inner(corr_op, psi0)

    println("Ground state energy: $e0")
    println("Ground state correlator: $corr")
    println("State vector length = ", size(psi0))
end

# use translation symmetry (less memory consuming)
function gs_correlator_translation_sym(N::Int, H::OpSum, corr_op::OpSum)
 
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
        sym_block = Spinhalf(N, irreps[k]) 
        lanczos_result = eigs_lanczos(H, sym_block)
        gs = lanczos_result.eigenvalues[1]
        if  gs < e0
            e0 = gs
            psi0 = lanczos_result.eigenvectors
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




