using LinearAlgebra
using Plots
using Combinatorics
# using BenchmarkTools
using Kronecker
using LaTeXStrings
using Arpack
# using KernelDensity
using Interpolations
using SparseArrays
# using ArnoldiMethod
using KrylovKit
using JLD2
using HDF5
using Printf
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2, 
    framestyle=:box, 
    # xtickfont=font(18),
    label=nothing,
    left_margin=4Plots.mm,
    bottom_margin=2Plots.mm
    # ytickfont=font(18),
    # legendfont=font(18)
)



n_sites=32
n_seeds = 1
ks=["Gamma.C4.A", "Gamma.C4.B", "Gamma.C4.Ea", "Gamma.C4.Eb", "M.C4.A", "M.C4.B", "M.C4.Ea", "M.C4.Eb", "X.C2.A", "X.C2.B"]#, "Delta.C1.A", "Sigma.C1.A", "Z0.C1.A", "Z1.C1.A"]
seeds = [i for i=1:n_seeds]
nup_start=10
n_ups = [i for i=nup_start:(div(n_sites,2))]

J1=1.00
n_eigs = 10


for seed in seeds
    plot()
    mineig = 0
    eigvs = []
    for n_up=n_ups
        for k in ks
            f = h5open(@sprintf("./outfiles/seed.%d/outfile.square.%d.J1.%.2f.nup.%d.k.%s.seed.%d.h5",seed,n_sites,J1,n_up,k,seed), "r")
            eig = read(f["Eigenvalues"])[1:n_eigs]
            eigvs = append!(eigvs,eig)
            close(f)
        end
    end
    mineig=minimum(eigvs)
    eigvs = eigvs-mineig*ones(length(eigvs))
    energies2 = round.(eigvs, digits=6)
    E0 = unique(energies2)
    sort!(E0)
    Stot = similar(E0)
    for i=1:length(E0)
        mask = findall(x -> x == E0[i], energies2)
        Sz = abs.(div.(mask,length(ks)*n_eigs).+(nup_start) .- n_sites / 2)
        vals = Sz .* (Sz .+ 1)
        max_val, arg = findmax(vals)
        Stot[i] =abs(max_val)
    end
    plot(Stot,E0,seriestype=:scatter,mc = :blue,legend=false,xlabel=L"S_\mathrm{tot}(S_{tot}+1)",ylabel=L"E/J_1")
    plot!(Stot,0.15*Stot,color = :black,legend=false,ylims=(0,10))
    savefig(@sprintf("outfile.square.%d.J1.%.2f.seed.%d-n.pdf",n_sites,J1,seed))
end
