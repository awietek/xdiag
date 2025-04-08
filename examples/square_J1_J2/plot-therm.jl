using LinearAlgebra
using Plots
using Combinatorics
# using BenchmarkTools
using Kronecker
using LaTeXStrings
# using Arpack
# using KernelDensity
using Interpolations
using SparseArrays
# using ArnoldiMethod
# using KrylovKit
using JLD2
using HDF5
using Printf

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2, 
    framestyle=:box, 
    # xtickfont=font(18),
    label=nothing
    # ytickfont=font(18),
    # legendfont=font(18)
)



n_sites=32
n_seeds = 1
ks=["Gamma.C4.A", "Gamma.C4.B", "Gamma.C4.Ea", "Gamma.C4.Eb", "M.C4.A", "M.C4.B", "M.C4.Ea", "M.C4.Eb", "X.C2.A", "X.C2.B", "Delta.C1.A", "None0.C1.A", "None1.C1.A", "Sigma0.C1.A", "Sigma1.C1.A", "Sigma2.C1.A", "Z0.C1.A"]
kmult=[1,1,1,1,1,1,1,1,2,2,4,4,4,4,4,4,4]
seeds = [i for i=2:2]
n_ups = [i for i=0:div(n_sites,2)]

J1=1.00
J2s=[i/10. for i=0:10]
Jchi = 0.0

temps = collect(range(0.01, 100, length=400))#2.5

mkpath("./data")
mkpath("./plot-therm")
mkpath("./plot-therm/spec-heat/")
mkpath("./plot-therm/susceptibility/")
mkpath("./plot-therm/entropy/")
mkpath("./plot-therm/partition/")


for J2 = J2s
    plot(title=L"J_2/J_1 = "*string(J2))
    SpecHeat = zeros(n_seeds,length(temps))
    partitions = zeros(n_seeds,length(temps))
    Susceptibility = zeros(n_seeds,length(temps))
    Energy = zeros(n_seeds,length(temps))
    Entropy = zeros(n_seeds,length(temps))
    n_seed=0
    for seed = seeds
        n_seed+=1
        f = h5open(@sprintf("outfile.square.%d.J1.%.2f.J2.%.2f.Jchi.%.2f.nup.%d.k.%s.seed.%d.h5",n_sites,J1,J2,Jchi,div(n_sites,2),"Gamma.C4.A",seed), "r")
        eig0 = read(f["Eigenvalues"])[1] # this is the ground state energy
        for ntemp = 1:length(temps)
            temp = temps[ntemp]
            partition = 0.0
            energy = 0.0
            energy2 =0.0
            sz2 =0.0
            # entropy = 0.0
            for nk = 1:length(ks)
                for n_up = n_ups
                    factor =2.0
                    if n_up == div(n_sites,2)
                        factor=1.0
                    end
                    f = h5open(@sprintf("/home/ssarkar/research/ChiralSpinHall/MyThirdCode/outfiles/seed.%d/outfile.square.%d.J1.%.2f.J2.%.2f.Jchi.%.2f.nup.%d.k.%s.seed.%d.h5",seed,n_sites,J1,J2,Jchi,n_up,ks[nk],seed), "r")
                    alphas = read(f["Alphas"])
                    betas = read(f["Betas"])
                    dims = read(f["Dimension"])
                    tmat = SymTridiagonal(alphas, betas[1:length(alphas)-1])
                    if length(alphas)==1
                        partition += exp(-(alphas[1]-eig0)/temp)*kmult[nk]*factor
                        energy += alphas[1]*exp(-(alphas[1]-eig0)/temp)*kmult[nk]*factor
                        energy2 += alphas[1]^2*exp(-(alphas[1]-eig0)/temp)*kmult[nk]*factor
                    elseif length(alphas)>1
                        F = eigen(tmat)
                        eig = F.values
                        vecs = F.vectors
                        partition += dims*sum([exp(-(eig[m]-eig0)/temp)*vecs[1,m]^2 for m=1:length(eig)])*kmult[nk]*factor
                        energy += dims*sum([eig[m]*exp(-(eig[m]-eig0)/temp)*vecs[1,m]^2 for m=1:length(eig)])*kmult[nk]*factor
                        energy2 += dims*sum([eig[m]^2*exp(-(eig[m]-eig0)/temp)*vecs[1,m]^2 for m=1:length(eig)])*kmult[nk]*factor
                        sz2 += dims*sum([(n_up-n_sites/2)^2*exp(-(eig[m]-eig0)/temp)*vecs[1,m]^2 for m=1:length(eig)])*kmult[nk]*factor
                    end
                end
            end
            SpecHeat[n_seed,ntemp] = (energy2/partition-(energy/partition)^2)/temp^2
            partitions[n_seed,ntemp] = partition
            Susceptibility[n_seed,ntemp] = sz2/partition/temp/n_sites
            Energy[n_seed,ntemp] = (energy/partition-eig0)
            Entropy[n_seed,ntemp] = (log(partition)+(energy/partition-eig0)/temp)/n_sites
        end
        plot!(temps,SpecHeat[n_seed,:],left_margin=15Plots.mm,bottom_margin=10Plots.mm,xlabel=L"T/J_1",ylabel=L"C_v")
    end
    savefig(@sprintf("plot-therm/spec-heat/outfile.square.%d.J1.%.2f.J2.%.2f.png",n_sites,J1,J2))
    plot(title=L"J_2/J_1 = "*string(J2))
    n_seed=0
    for seed = seeds
        n_seed+=1
        plot!(temps,partitions[n_seed,:],left_margin=20Plots.mm,bottom_margin=10Plots.mm,xlabel=L"T/J_1",ylabel=L"Z")
    end
    savefig(@sprintf("plot-therm/partition/outfile.square.%d.J1.%.2f.J2.%.2f.png",n_sites,J1,J2))
    plot(title=L"J_2/J_1 = "*string(J2))
    n_seed=0
    for seed = seeds
        n_seed+=1
        plot!(temps[1:200],Susceptibility[n_seed,1:200], xlims=(1e-2, temps[200]), ylims=(1e-5, 1.2e-1),left_margin=20Plots.mm,bottom_margin=10Plots.mm,xlabel=L"T/J_1",ylabel=L"\chi")
    end
    savefig(@sprintf("plot-therm/susceptibility/outfile.square.%d.J1.%.2f.J2.%.2f..png",n_sites,J1,J2))
    plot(title=L"J_2/J_1 = "*string(J2))
    n_seed=0
    for seed = seeds
        n_seed+=1
        plot!(temps[1:200],Entropy[n_seed,1:200], xlims=(1e-2, temps[200]),left_margin=20Plots.mm,bottom_margin=10Plots.mm,xlabel=L"T/J_1",ylabel=L"s")
    end
    savefig(@sprintf("plot-therm/entropy/outfile.square.%d.J1.%.2f.J2.%.2f.png",n_sites,J1,J2))
    fid = h5open(@sprintf("data/outfile.square.%d.J1.%.2f.J2.%.2f.h5",n_sites,J1,J2),"w")
    fid["SpecHeat"] = SpecHeat
    fid["partitions"] = partitions
    fid["Susceptibility"] = Susceptibility
    fid["Energy"] = Energy
    fid["Entropy"] = Entropy
    close(fid)
end



plot()
let n_seed = 0
    for seed = seeds
        n_seed+=1
        for J2 = [i/10. for i=[1,3,5,7,9,10]]
            fid = h5open(@sprintf("data/outfile.square.%d.J1.%.2f.J2.%.2f.h5",n_sites,J1,J2),"r")
            entropy = read(fid["Entropy"]);
            susceptibility = read(fid["Susceptibility"])
            partitions=read(fid["partitions"])
            id=1
            j=1
            while partitions[n_seed,j]<10
                j=j+1
            end
            id=j
            wlratio = 4*pi^2*(temps.*susceptibility[n_seed,:])./(3*entropy[n_seed,:])
            plot!(temps[id:200],wlratio[id:200], xlims=(1e-2, temps[200]),ylims=(0.0,3.0),left_margin=20Plots.mm,bottom_margin=10Plots.mm,xlabel=L"T",ylabel=L"\mathrm{Wilson\ ratio\ }R",label=L"J_2 = "*string(J2), dpi=1000)
        end
        savefig(@sprintf("plot-therm/WR.square.%d.J1.%.2f.Jchi.%.2f.png",n_sites,J1,Jchi))
    end
end
