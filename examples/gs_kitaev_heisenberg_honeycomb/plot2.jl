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
    left_margin=8Plots.mm,
    bottom_margin=2Plots.mm
    # ytickfont=font(18),
    # legendfont=font(18)
)
colors = palette(:default)#palette(:thermal,length(14:div(n_sites,2))+1)
markers = filter((m->begin
                m in Plots.supported_markers()
            end), Plots._shape_keys)


n_sites=24
n_seeds = 2
ks=["Gamma.C2.A", "Gamma.C2.B", "M0.C2.A", "M0.C2.B", "M1.C2.A", "M1.C2.B", "M2.C2.A", "M2.C2.B", "K0.C1.A", "K1.C1.A", "X0.C1.A", "X1.C1.A", "X2.C1.A"]
ksl=["Γ.C2.A", "Γ.C2.B", "M0.C2.A", "M0.C2.B", "M1.C2.A", "M1.C2.B", "M2.C2.A", "M2.C2.B", "K0.C1.A", "K1.C1.A", "X0.C1.A", "X1.C1.A", "X2.C1.A"]
seeds = [i for i=2:n_seeds]
Ks = [i/10. for i=-40:0]

J=1.00

for seed in seeds
    plot(xlabel=L"K/J",ylabel=L"E/J",ylims=(-0.02,2.0),xlims=(-0.02,5.0),xticks=(0.0:0.4:4,[@sprintf("%.2f",i) for i=0.0:0.4:4]),dpi=400)
    for K=Ks
        mineig = 100.
        for k in ks
            f = h5open(@sprintf("outfiles/seed.%d/outfile.honeycomb.%d.J.%.2f.KX.%.2f.KY.%.2f.KZ.%.2f.k.%s.seed.%d.h5",seed,n_sites,J,K,K,K,k,seed), "r")
            eigs =read(f["Eigenvalues"])
            if mineig>eigs[1]
                mineig = eigs[1]
            end
        end
        c=0
        for k in ks
            c=c+1
            f = h5open(@sprintf("outfiles/seed.%d/outfile.honeycomb.%d.J.%.2f.KX.%.2f.KY.%.2f.KZ.%.2f.k.%s.seed.%d.h5",seed,n_sites,J,K,K,K,k,seed), "r")#,J,K,K,K
            eig = read(f["Eigenvalues"]).-mineig
            KK=abs(K)*[1 for i=1:length(eig)]
            if K==Ks[1]
                plot!(KK,eig,seriestype=:scatter,m=markers[c],mc=colors[c],label=ksl[c])
            else
                plot!(KK,eig,seriestype=:scatter,m=markers[c],mc=colors[c],primary=false)
            end
            close(f)
        end
    end
    savefig(@sprintf("outfile.kitaev.%d.seed.%d-n.png",n_sites,seed))
end
