using XDiag
using LinearAlgebra
using HDF5
using Plots
using IterTools

# pythonplot()

function poles_weights(outfile::String; M::Int64=-1)
    e0 = h5read(outfile, "e0")
    nrm, alphas, betas, eigs = h5read.((outfile,), ["norm", "alphas", "betas", "eigs"])
    if M > 0
        try
            alphas = alphas[1:M]
            betas = betas[1:M]
        catch e
            println(e)
        end
    end
    tmat = SymTridiagonal(alphas, betas)
    es, evecs = eigen(tmat)
    # @show es[1], eigs[1]

    poles = es .- e0
    weights = evecs[1, :] .^ 2 .* nrm^2
    sortdisp(poles, weights/nrm^2)
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

function sortdisp(p,w;max=10)
    upto=min(max,length(p))
    display([p[sortperm(w,rev=true)][1:upto] sort(w,rev=true)[1:upto]])
end

let
    source = "../../misc/data/examples_output/tJ_conductivity_jll.h5"
    omegas = collect(range(0., 10.0, length=40000))
    # eta = omegas[2]-omegas[1]
    eta = 0.1
    fig = plot()
    for M in [30, 50, 100, 200]
        p, w = poles_weights(source; M=M)
        s = spectrum(p, w, omegas, eta)
        plot!(omegas, s, label="M=$M")
    end
    savefig(fig, "cond_jll.pdf")
end
