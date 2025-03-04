using XDiag
using LinearAlgebra
using HDF5
using Plots

function plot_sf(outfile::String, omegas::Vector{Float64}, eta)
    alls = Matrix{Float64}(undef, 16, length(omegas))
    for k in 0:15
        p, w = poles_weights(outfile, k)
        s = spectrum(p, w, omegas, eta)
        alls[k+1, :] = s
    end
    return heatmap(collect(range(0, 15, step=1)), omegas, alls',
        xlabel="\$k\$", ylabel="\$\\omega\$",
        interpolation=:false)
end

function poles_weights(outfile::String, k::Int64)
    e0 = h5read(outfile, "/e0")
    nrm, alphas, betas, eigs = h5read.((outfile,), ["/$k/norm", "/$k/alphas", "/$k/betas", "/$k/eigs"])
    tmat = SymTridiagonal(alphas, betas[1:end-1])
    es, evecs = eigen(tmat)
    # es == eigs

    poles = es .- e0
    weights = evecs[1, :] .^ 2 .* nrm^2
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

let
    source = "../../misc/data/examples_output/spinhalf_chain_structure_factor.h5"
    omegas = collect(range(0.0, 3.0, length=1000))
    eta = 10 * (omegas[2] - omegas[1])
    fig = plot_sf(source, omegas, eta)
    savefig(fig, "spekt.pdf")
end
