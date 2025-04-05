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


n_sites=36
n_seeds = 1
ks=["Gamma.C6.A", "Gamma.C6.B", "Gamma.C6.E1a", "Gamma.C6.E1b", "Gamma.C6.E2a", "Gamma.C6.E2b", "K.C3.A", "K.C3.Ea", "K.C3.Eb", "M.C2.A", "M.C2.B"]#, "Z.C1.A", "Delta.C1.A", "Sigma.C1.A", "Z0.C1.A", "Z1.C1.A"]
ksl=["Γ.C6.A", "Γ.C6.B", "Γ.C6.E1a", "Γ.C6.E1b", "Γ.C6.E2a", "Γ.C6.E2b", "K.C3.A", "K.C3.Ea", "K.C3.Eb", "M.C2.A", "M.C2.B"]
seeds = [i for i=1:n_seeds]
nup_start=10
n_ups = [i for i=nup_start:(div(n_sites,2))]#

J1=1.00
J2=-1.0
J3=-2.0
n_eigs = 10

for seed in seeds
    plot(xlabel=L"S_z",ylabel=L"E/J_1",ylims=(-0.02,1),xlims=(-0.1,10),dpi=400,xticks=(0:2:10,[string(i) for i=0:2:10]))
    mineig = 0
    eigvs = []
    Sz = []
    f = h5open(@sprintf("outfiles/seed.%d/outfile.kagome.%d.J1.%.2f.J2.%.2f.J3.%.2f.nup.%d.k.%s.seed.%d.h5",seed,n_sites,J1,J2,J3,div(n_sites,2),"Gamma.C6.A",seed), "r")
    mineig = read(f["Eigenvalues"])[1]
    for n_up=n_ups
        c=0
        for k in ks
            c=c+1
            f = h5open(@sprintf("outfiles/seed.%d/outfile.kagome.%d.J1.%.2f.J2.%.2f.J3.%.2f.nup.%d.k.%s.seed.%d.h5",seed,n_sites,J1,J2,J3,n_up,k,seed), "r")
            eig = read(f["Eigenvalues"])[1:n_eigs].-mineig
            sz=abs(n_up- n_sites / 2)*[1 for i=1:n_eigs]
            if n_up==nup_start
                plot!(sz,eig,seriestype=:scatter,m=markers[c],mc=colors[c],label=ksl[c])
            else
                plot!(sz,eig,seriestype=:scatter,m=markers[c],mc=colors[c],primary=false)
            end
            close(f)
        end
    end
    savefig(@sprintf("outfile.kagome.%d.J1.%.2f.J2.%.2f.J3.%.2f.seed.%d-n.png",n_sites,J1,J2,J3,seed))
end
