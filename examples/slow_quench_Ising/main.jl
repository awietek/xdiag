using XDiag
using LinearAlgebra
using Printf
using GLMakie
using LaTeXStrings

# Build a LaTeX legend label for τ_Q (used in plotting)
function latex_label(tauQ::Union{Float64, Int})
    return LaTeXString("\\tau_Q = " * string(tauQ))
end

# Construct the Ising Hamiltonian with transverse and longitudinal fields
function make_hamiltonian(N::Int, hx::Float64, hz::Float64)
    op_sum = OpSum()

    # Nearest-neighbor Ising term: -σᶻᵢ σᶻᵢ₊₁
    for i in 1:(N-1)
        key_zz = "zz_$i"
        op_sum += key_zz * Op("SzSz", [i, i+1])
        op_sum[key_zz] = -4.0  # compensate for Sz = ±1/2
    end

    # Transverse field: -hx * σˣᵢ = -hx (S⁺ + S⁻)
    for i in 1:N
        op_sum += "sp_$i" * Op("S+", [i])
        op_sum["sp_$i"] = -hx
        op_sum += "sm_$i" * Op("S-", [i])
        op_sum["sm_$i"] = -hx
    end

    # Longitudinal field: -hz * σᶻᵢ
    for i in 1:N
        op_sum += "hz_$i" * Op("Sz", [i])
        op_sum["hz_$i"] = -2 * hz
    end

    return op_sum
end

# Measure longitudinal magnetization m_z
function measure(psi::State, N::Int)
    mz = 0.0
    for i in 1:N
        mz += 2.0 * real(inner(Op("Sz", [i]), psi))  # convert Sz = ±1/2 to ±1
    end
    return mz / N
end

# Compute ground state for given parameters
function ground_state(N::Int, hx::Float64, hz::Float64)
    block = Spinhalf(N)
    H = make_hamiltonian(N, hx, hz)
    return eig0(H, block)
end

function main()
    # System and ramp parameters
    N = 8
    hx = 0.2
    hz_in, hz_fin = -3.0, 3.0
    dt = 0.05
    tauQ_list = [64, 256, 1024]  # different ramp speeds

    results = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()

    # Simulate slow longitudinal field ramps for various τ_Q
    for τQ in tauQ_list
        @printf("\n=== Running quench for tauQ=%.1f ===\n", τQ)
        total_time = (hz_fin - hz_in) * τQ
        nsteps = Int(ceil(total_time / dt))
        dt_eff = total_time / nsteps

        e0, psi0 = ground_state(N, hx, hz_in)

        hz_vals = zeros(nsteps+1)
        mz_vals = zeros(nsteps+1)
        hz_vals[1] = hz_in
        mz_vals[1] = measure(psi0, N)

        psi_t = psi0
        for step in 1:nsteps
            t_cur = step * dt_eff
            frac = t_cur / total_time
            hz_now = hz_in + (hz_fin - hz_in) * frac

            H = make_hamiltonian(N, hx, hz_now)
            psi_t = time_evolve(H, psi_t, dt_eff)

            hz_vals[step+1] = hz_now
            mz_vals[step+1] = measure(psi_t, N)

            if step % 500 == 0
                @printf("step=%5d/%5d => t=%.2f, h_z=%.2f, m_z=%.3f\n",
                        step, nsteps, t_cur, hz_now, mz_vals[step+1])
            end
        end
        results[τQ] = (hz_vals, mz_vals)
    end

    # Plotting m_z vs h_z(t) for each ramp rate
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1];
        title="Bubble Nucleation through Slow Quench",
        xlabel=L"$h_z(t)$", ylabel=L"$m_z$",
        xticks = ([-3.0, -2, -1,  0.0, 1, 2, 3.0],
                  [L"-3", L"-2", L"-1", L"0", L"1", L"2", L"3"]),
        
        yticks = ([-1.0, -0.5, 0.0, 0.5],
                  [L"-1", L"-0.5", L"0", L"0.5"]),
        xticksize=20, yticksize=20,
        xticklabelsize=36, yticklabelsize=36,
        xlabelsize=36, ylabelsize=36,
        titlesize=29
    )

    colors = (:blue, :red, :green)
    for (i, τQ) in enumerate(tauQ_list)
        hzvals, mzvals = results[τQ]
        lines!(ax, hzvals[1:5:end], mzvals[1:5:end];
            color=colors[i], linewidth=3,
            label=latex_label(τQ)
        )
    end

    axislegend(ax, position=:lt, labelsize=40)
    save("magnetization_slow_quench.png", fig)
    println("\n Plot saved to 'magnetization_slow_quench.png'")
end

main()
