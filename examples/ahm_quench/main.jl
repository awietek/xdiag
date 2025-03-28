using XDiag
using CairoMakie

"""Struct for lattice sites"""
struct LatticeBond
    s1::Int64
    x1::Float64
    y1::Float64
    s2::Int64
    x2::Float64
    y2::Float64
end

"""Returns list of all bonds on the square lattice"""
function square_lattice(Lx::Int,Ly::Int;xperiodic=true,yperiodic=true)

    bonds = Vector{LatticeBond}()

    for x1 in 1:Lx, y1 in 1:Ly
        s1 = (x1-1)*Ly+y1
        for (x2,y2) in [( mod(x1,Lx)+1, y1 ),( x1, mod(y1,Ly)+1 )]
            s2 = (x2-1)*Ly+y2
            set = reduce(hcat,[[s1,x1,y1],[s2,x2,y2]])
            set = set[:,sortperm(set[1,:])]     # sort according to s1 and s2
            b = LatticeBond(set[:,1]...,set[:,2]...)
            push!(bonds,b)
        end
    end
    return bonds
end

"""Returns the Hamiltonian for the Hubbard model"""
function ham(U::Float64,Lx::Int64,Ly::Int64)

    lattice = square_lattice(Lx,Ly)
    ops = OpSum()
    
    for b in lattice # hopping term
        ops += "t" * Op("Hop", [b.s1,b.s2])     # both for spin up and down
    end

    ops += "U" * Op("HubbardU")     # applies Hubbard interaction over the entire lattice
    for i in 1:Lx*Ly
        ops += "mu" * Op("Ntot",i)     # chemical potental
    end
    
    ops["t"] = 1.0
    ops["U"] = U
    #ops["mu"] = -U/2
    ops["mu"] = 0.0
    
    return ops

end

let
    U = 0.0
    Lx,Ly = (4,3)
    N = Lx*Ly
        # ground state for U
    ops = ham(U,Lx,Ly)
    hspace = Electron(N, N ÷ 2, N ÷ 2)

    e,psi = eig0(ops,hspace)

    U = -10.0
    ops = ham(U,Lx,Ly)
    dt = 0.1
    time = 10

    times = range(dt,time,length=Int(time/dt))
    obs_vec = Array{Float64}(undef,length(times))

    obs = OpSum()
    obs += Op("HubbardU")

    for i in 1:length(times)
        time_evolve_inplace(ops,psi,float(dt))
        # do measurements
        obs_vec[i] = real(inner(obs,psi))
    end

    f = Figure()
    ax = Axis(f[1,1],
        xlabel=L"time, $t$",
        ylabel=L"$\sum_i \langle n_{i\uparrow} n_{i\downarrow} \rangle $"
        )


    lines!(ax,times,obs_vec)
    save("../plots/quench/Lx($Lx)_Ly($Ly).png",f)
end
