using XDiag
using Plots 

function main()

    # specify system parameters
    N = 18
    J = 1.0

    # fix number of up spins
    Nup = N÷2

    # build Hamiltonian
    H = OpSum()
    for i in 1:N
        H += "J" * Op("SdotS", [i,mod1(i+1, N)])
    end
    H["J"] = J
    
    # set up translation symmetry
    T = Permutation(circshift(1:N, -1))
    C_N = PermutationGroup([T^k for k in 0:(N-1)])
    character_table = [ [C_N_character(N, k, p) for p in 0:(N-1)] for k in 0:(N-1)]
    irreps = [ Representation(C_N, character_table[k]) for k in 1:N ]   

    # sort states in Sz_tot = 0 sector by translation irrep
    lvl_per_block = 5
    TOS_x_vals = Vector{Float64}(undef, N*lvl_per_block)
    TOS_y_vals  = Vector{Float64}(undef, N*lvl_per_block)
    for k in eachindex(irreps)
        sym_block = Spinhalf(N, Nup, irreps[k])
        # set neigvals = ... for reasonable convergence of lowest eigenvalues
        lanczos_result = eigvals_lanczos(H, sym_block, neigvals=3) 
        for i in 1:lvl_per_block
            TOS_x_vals[(k-1)*lvl_per_block + i] = 2*pi*(k-1)/N
            TOS_y_vals[(k-1)*lvl_per_block + i] = lanczos_result.eigenvalues[i]
        end
    end

    # create TOS plot
    scatter(
        TOS_x_vals,
        TOS_y_vals,
        title="Tower of States, AFM Heisenberg Chain, N=$N\n",
        xlabel="lattice momentum k",
        ylabel="lowest eigenstates",
        legend=false,
        xticks = ([0:π/2:2*π;], ["0", "\\pi/2", "\\pi", "3\\pi/2", "2\\pi"]),
        right_margin=5Plots.mm,
        dpi=300)
end


function C_N_character(N::Int, k::Int, p::Int)
    return exp( im * 2 * pi * p * k * 1.0 / N )
end


main()
