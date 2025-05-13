using XDiag
using LinearAlgebra
using HDF5
using Plots
using PyPlot
using IterTools

# pythonplot()

function plot_dos(outfile::String, omegas::Vector{Float64}, eta::Float64=0.15)
    f = h5open(outfile)
    @show ikxs = filter(s -> occursin(r"^\d+$", s), keys(f)) # assuming ikxs == ikys
    alls = zeros(Float64, length(omegas))
    sumRule = 0.0
    for (ikx, iky) in IterTools.product(ikxs, ikxs)
        sumRule += h5read(outfile, "/$ikx/$iky/norm")^2 # should divide by nsites
        p, w = poles_weights(outfile, "$ikx/$iky")
        s = spectrum(p, w, omegas, eta)
        alls += s
    end
    display("sum rule: $sumRule")

    return plot(omegas, alls, xlabel="\$\\omega\$", ylabel="DOS")
end

function plot_sf(outfile::String, omegas::Vector{Float64}, eta::Float64=0.15)
    irreps = Dict{String,Vector{Float64}}(
        "Gamma.C1.A" => [0.0, 0.0],
        "M.C1.A" => [3.1415926535897931, 3.1415926535897931],
        "Sigma0.C1.A" => [1.5707963267948966, 1.5707963267948966],
        "Sigma1.C1.A" => [1.5707963267948966, -1.5707963267948966],
        "Sigma2.C1.A" => [-1.5707963267948966, 1.5707963267948966],
        "Sigma3.C1.A" => [-1.5707963267948966, -1.5707963267948966],
        "X0.C1.A" => [3.1415926535897931, 0.0000000000000000],
        "X1.C1.A" => [0.0000000000000000, 3.1415926535897931])
    alls = zeros(Float64, length(irreps), length(omegas))
    alls2 = zeros(Float64, length(irreps), length(omegas))
    sumRule = 0.0
    for (k, (name, momentum)) in enumerate(irreps)
        sumRule += h5read(outfile, "/$name/irrep/norm")^2 # should divide by nsites
        p, w = poles_weights(outfile, name, "irrep")
        s = spectrum(p, w, omegas, eta)
        alls[k, :] = s

        p, w = poles_weights(outfile, name, "tbc")
        s = spectrum(p, w, omegas, eta)
        alls2[k, :] = s
        display(plot(omegas, hcat(alls[k,:], alls2[k,:])))
        sleep(2)
    end
    display("sum rule: $sumRule")
    @show size(alls), length(omegas), length(irreps),
    # return heatmap(range(0, length(irreps) - 1, step=1), omegas, alls',
    #     xlabel="\$k\$", ylabel="\$\\omega\$")
    return plot(omegas, hcat(alls',alls2'), xlabel="\$\\omega\$")
end

function plot_momentum_cut(outfile::String, omegas::Vector{Float64}, eta::Float64=0.15)
    f = h5open(outfile)
    @show ikxs = filter(s -> occursin(r"^\d+$", s), keys(f)) # assuming ikxs == ikys
    alls = zeros(Float64, length(ikxs), length(omegas))
    sumRule = 0.0
    for (k, ikx) in enumerate(ikxs)
        for component in ["ccdag", "cdagc"]
            sumRule += h5read(outfile, "/$ikx/$ikx/$component/norm")^2 # should divide by nsites
            p, w = poles_weights(outfile, "$ikx/$ikx", component)
            s = spectrum(p, w, omegas, eta)
            alls[parse(Int, ikx)+1, :] += s
        end
    end
    display("sum rule: $sumRule")
    display("$(size(alls)), $(length(omegas)), $(length(ikxs))")

    # display(alls)
    fig = heatmap(range(0, pi, length=length(ikxs)), omegas, (alls'),
        # clim=(-5,1),
        xlabel="\$k\$", ylabel="\$\\omega\$", xticks=([0, pi/4, pi/2, 3*pi/4, pi], ["0","π/4","π/2","3π/4","π"]),
        framestyle=:box)
    # return plot(omegas, alls')
    ax = PyPlot.gca()
    ax.tick_params(top=true, right=true)
    return fig
end

function poles_weights(outfile::String, name::String)
    e0 = h5read(outfile, "/e0")
    nrm, alphas, betas, eigs = h5read.((outfile,), ["/$name/norm", "/$name/alphas", "/$name/betas", "/$name/eigs"])
    @show name, nrm
    tmat = SymTridiagonal(alphas, betas)
    es, evecs = eigen(tmat)
    # @show es[1], eigs[1]

    poles = e0 .- es
    weights = evecs[1, :] .^ 2 .* nrm^2
    sortdisp(poles, weights / nrm^2)
    return poles, weights
end

function poles_weights(outfile::String, name::String, component::String)
    e0 = h5read(outfile, "/e0")
    nrm, alphas, betas, eigs = h5read.((outfile,), ["/$name/$component/norm", "/$name/$component/alphas", "/$name/$component/betas", "/$name/$component/eigs"])
    @show name, nrm
    tmat = SymTridiagonal(alphas, betas)
    es, evecs = eigen(tmat)
    # @show es[1], eigs[1]

    poles = e0 .- es
    weights = evecs[1, :] .^ 2 .* nrm^2
    sortdisp(poles, weights / nrm^2)
    return poles, weights
end

function spectrum(poles::Vector{Float64}, weights::Vector{Float64},
    omegas::Vector{Float64}, eta::Float64=0.1; cutoff::Float64=0.0, filter::Float64=1e-6)
    fpoles = [poles[i] for i in 1:length(weights) if abs(weights[i]) > filter]
    fweights = [weights[i] for i in 1:length(weights) if abs(weights[i]) > filter]
    diffs = omegas .- fpoles'
    gaussians = exp.(-(diffs ./ (2.0 * eta)) .^ 2.0) ./ (eta * sqrt(2.0 * pi))
    if abs(cutoff) > 1e-12
        idxs = findall(abs.(fpoles) .> cutoff)
        s = zeros(length(omegas))
        for idx in idxs
            s += gaussians[:, idx] * fweights[idx]
        end
        return s
    else
        return gaussians * fweights
    end
end

function sortdisp(p, w; max=10)
    upto = min(max, length(p))
    display([p[sortperm(w, rev=true)][1:upto] sort(w, rev=true)[1:upto]])
end

let
    source = "../../misc/data/examples_output/hubbard_greens_f.h5"
    # omegas = collect(range(-8, 8, length=4000))
    # eta = 5 * (omegas[2] - omegas[1])
    pyplot()

    # fig = plot_dos(source, omegas, 0.05)
    # savefig(fig, "dos.pdf")

    # fig = plot_sf(source, omegas, 0.15)
    # display(fig)
    # sleep(100)
    # savefig(fig, "spekt.pdf")

    omegas = collect(range(-4.0, 6.0, length=2000))
    fig = plot_momentum_cut(source, omegas)
    # display(fig)
    # sleep(10)
    Plots.savefig(fig, "spekt_momentum_cut.pdf")
end
