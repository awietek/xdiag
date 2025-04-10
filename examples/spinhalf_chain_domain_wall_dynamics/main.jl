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
    psi0_vec = vcat(repeat(["Up"], N÷2), repeat(["Dn"], N÷2))
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

    # plot Sz expectation value
    @show heatmap(
        Sz_expectation,
        c = :thermal,
        xlabel="chain coordinate",
        ylabel="time step",
        title = "Sz expectation value over time")
    
end


main()
