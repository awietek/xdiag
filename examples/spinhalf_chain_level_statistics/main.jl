using XDiag
using LinearAlgebra # exact diagonalization
using Plots # optional: plot histograms


function main()

    N = 18

    # definition of integrable model
    H_i = OpSum()
    for i in 1:N
        H_i += "J" * Op("SzSz", [i,mod1(i+1, N)])
        H_i += "Delta" * Op("Exchange", [i, mod1(i+1, N)])
    end

    # definition of non-integrable model
    H_ni = OpSum()
    for i in 1:N
        H_ni += "J" * Op("SzSz", [i, mod1(i+1, N)])
        H_ni += "Delta" * Op("Exchange", [i, mod1(i+1, N)])
        H_ni += "J2" * Op("SzSz", [i, mod1(i+2, N)])
    end

    # assign coupling values
    J = 1.0
    Delta = J2 = 0.5

    H_i["J"] = H_ni["J"] = J
    H_i["Delta"] = H_ni["Delta"] = Delta
    H_ni["J2"] = J2

    # compute level statistics (remember to eliminate trivial symmetries!)
    H_i_statistics = compute_level_statistics(N, H_i)
    H_ni_statistics = compute_level_statistics(N, H_ni)

    # optional plot of histograms (only in Julia version!)
    plot_histograms(H_i_statistics, H_ni_statistics)
end


function compute_level_statistics(N::Int, H::OpSum) :: Vector{Float64}

    # fix magnetization (Sz_tot = 0 still has spin-flip symmetry!)
    Nup = N÷2 + 1

    # fix lattice momentum (k = 0, N/2 still have parity symmetry!)
    k = 1
    T = Permutation(circshift(1:N, -1))
    C_N = PermutationGroup([T^p for p in 0:(N-1)])
    irrep_k_characters = [C_N_character(N, k, p) for p in 0:(N-1)]
    irrep_k = Representation(C_N, irrep_k_characters)

    # block of Hamiltonian without remaining symmetries
    block = Spinhalf(N, Nup, irrep_k) 

    # find its eigenspectrum
    Hmat = matrix(H, block)
    eigenvalues = eigvals(Hermitian(Hmat))

    # compute level statistics (taking only inner most half of spectrum)
    N_levels = size(block)
    s_start = trunc(Int, 0.25*N_levels)
    s_stop = trunc(Int, 0.75*N_levels)
    s_num = s_stop - s_start
    s_arr = Vector{Float64}(undef, s_num)
    s_arr_sum = 0.0
    for i in 1:s_num
        s_arr[i] = eigenvalues[s_start+i+1] - eigenvalues[s_start+i]  
        s_arr_sum += s_arr[i]
    end

    # normalize
    return s_arr / (s_arr_sum / s_num)
end 


function C_N_character(N::Int, k::Int, p::Int)
    return exp( im * 2 * pi * p * k * 1.0 / N )
end


function plot_histograms(integ_stat::Vector{Float64}, noninteg_stat::Vector{Float64})
    # plot histograms
    smax = 3
    Nbins = 20
    bins = LinRange(0, smax, Nbins)
    histogram(
        [integ_stat noninteg_stat],
        bins=bins,
        normalize=:pdf,
        fillalpha=0.3,
        label=["integrable" "non-integrable"],
        xlabel="s",
        ylabel="P(s)"
    )

    # plot Wigner Dyson distribution
    x_vals = LinRange(0, smax, 200)
    WD_y_vals = WD_func.(x_vals)
    plot!(x_vals, WD_y_vals, color=:red, linewidth=3, label="Wigner-Dyson")

    # plot Poisson distribution
    Pois_y_vals = Poisson_func.(x_vals)
    plot!(x_vals, Pois_y_vals, color=:black, linewidth=3, label="Poisson")
end


function WD_func(s::Float64) :: Float64
    return (pi*s/2) * exp(-pi*s^2/4)
end


function Poisson_func(s::Float64) :: Float64
    return exp(-s)
end




main()