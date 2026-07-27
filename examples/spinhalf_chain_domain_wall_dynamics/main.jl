using XDiag
using Plots 

function main()

    # define open ferromagnetic XXZ chain
    N = 16
    J = 0.1
    Delta = 0.5

    H = OpSum()
    for i in 1:(N-1)
        H += "J" * Op("SzSz", [i, i+1])
        H += "Delta" * Op("Exchange", [i, i+1])
    end

    H["J"] = J
    H["Delta"] = Delta

    # define initial state with domain wall
    block = Spinhalf(N)
    psi0_vec = vcat(repeat([1], N÷2), repeat([0], N÷2))
    psi = product_state(block, psi0_vec)

    # time evolve and measure Sz expectation value
    dt = 0.5
    Nt = 30
    Sz_expectation = Matrix{Float64}(undef, Nt, N)
    for t_step in 1:Nt
        time_evolve_inplace(H, psi, dt)
        for i in 1:N
            Sz_expectation[t_step, i] = real(inner(Op("Sz", [i]), psi))
        end
    end

    # have some dummy output here to check against the C++ version of the example
    println("Computation of Sz expectation value over time successful!")
    println("First 10 entries of Sz_expectation at the final time step:")
    for i in 1:10
        println("$(i-1): ", Sz_expectation[Nt, i])
    end

    # plot Sz expectation value
    @show heatmap(
        Sz_expectation,
        c = :thermal,
        xlabel="chain coordinate",
        ylabel="time step",
        title = "Sz expectation value over time")
    
    return 0;
end


main()


#=
Expected output: ----------

Computation of Sz expectation value over time successful!
First 10 entries of Sz_expectation at the final time step:
0: 0.4432301014161557
1: 0.3949290171366142
2: 0.26315340423182687
3: 0.19962207952811017
4: 0.19857671600477683
5: 0.1273102824169592
6: 0.06956175013570165
7: 0.039150756300690974
8: -0.039150756300690905
9: -0.06956175013570176
=#

