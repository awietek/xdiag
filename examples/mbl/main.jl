using XDiag
using CairoMakie
using LinearAlgebra
using Statistics

function diagonalise(L,W,J,h)

    # construct Hilbert space
    hspace = Spinhalf(L)

    # construct Hamiltonian
    ops = OpSum()

    for i in 1:L-1
        ops += J[i] * Op("SzSz",[i,i+1])
    end

    for i in 1:L
        ops += h[i] * Op("Sz",i)
        ops += 0.5 * Op("S+",i)
        ops += 0.5 * Op("S-",i)
    end

    # Lanczos
    # res = eigs_lanczos(ops,hspace)
    H = matrix(ops, hspace)
    eig = eigen(Hermitian(H))

    return eig.values, eig.vectors

end

function level_spacing_ratio(eigs)

    diff = []
    for i in 1:length(eigs)-1
        push!(diff,eigs[i+1] - eigs[i])
    end

    r = []
    for i in 1:length(diff)-1
        push!(r,min(diff[i],diff[i+1])/max(diff[i],diff[i+1]))
    end

    return mean(r)

end

function ipr_func(state)

    ipr = 0
    for i in 1:length(state)
        ipr += abs(state[i])^4
    end

    return ipr

end

function main()
    # set model parameters
    L = 10
    Ws = range(0,step=0.05,length=100)
    nres = 20
    rs = Array{Float64}(undef,nres,length(Ws))
    f = Figure()
    
    neigs = 4   # number of low lying eigenstates to check IPR for
    iprs = zeros(nres,neigs,length(Ws))
    
    for ii in 1:nres
        J = (rand(Float64,L-1) .* 0.4) .+ 0.8   # J between 0.8 and 1.2

        for (jj,W) in enumerate(Ws)
            println("Now computing W=$W")
            h = (rand(Float64,L) .* 2*W) .- W   # h between -W and W

            evals, evecs = diagonalise(L,W,J,h)
            rs[ii,jj] = level_spacing_ratio(evals)
        
            for k in 1:neigs
                evec = evecs[:,k]
                iprs[ii,k,jj] = ipr_func(evec)
            end

        end

    end

    r = dropdims(mean(rs,dims=1),dims=1)
    ipr = dropdims(mean(iprs,dims=1),dims=1)

    ax1 = Axis(f[1,1],
        xlabel=L"$W$",
        ylabel="IPR",
        )

    for k in 1:neigs
        scatter!(ax1,Ws,ipr[k,:])
    end

    ax2 = Axis(f[1,2],
        xlabel=L"$W$",
        ylabel=L"$\langle r \rangle$",
        )

    scatter!(ax2,Ws,r)

    hlines!(ax2,0.386,label="Poisson",color=:red)
    hlines!(ax2,0.530,label="Gaussian",color=:green)

    axislegend(ax2)

    save("../plots/mbl/L($L).png",f)

end

main()