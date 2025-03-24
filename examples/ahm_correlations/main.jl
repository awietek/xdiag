
using XDiag
using CairoMakie
using HDF5

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
    
    ops["t"] = 1.0
    ops["U"] = U
    
    return ops

end

function main(U::Float64,Lx::Int64,Ly::Int64)
    
    say_hello()
    N = Lx*Ly   # number of particles, half-filling

    ops = ham(U,Lx,Ly) # create the Hamiltonian

    # Hilbert space with nup = ndn
    block = Electron(N, N ÷ 2, N ÷ 2)

    set_verbosity(0)
    e0, psi0 = eig0(ops,block)

    corr = zeros(N)
    for j in 2:N
        op = Op("NtotNtot",[1,j])
        corr[j] = inner(op,psi0)
    end

    corr = reshape(corr,(Lx,Ly))

    filename = "data/ahm_correlations/U($U)_Lx($Lx)_Ly($Ly).h5"

    h5open(filename,"w") do f
        write(f,"correlator",corr)
    end

    f = Figure()
    ax = Axis(f[1,1],
        xlabel="x",
        ylabel="y"
        )
    
    heatmap!(ax,corr)

    Colorbar(f[1, 2], limits = extrema([corr...]), colormap = :viridis,
    flipaxis = false)

    save("plots/ahm_correlations/U($U)_Lx($Lx)_Ly($Ly).png",f)

end

U = -10.0
Lx, Ly = (4,4)

main(U,Lx,Ly)