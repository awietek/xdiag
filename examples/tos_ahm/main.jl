using XDiag
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
    for i in 1:Lx*Ly
        ops += "mu" * Op("Ntot",i)     # chemical potential
    end
    
    ops["t"] = 1.0
    ops["U"] = U
    ops["mu"] = -U/2
    
    return ops

end

function main(U::Float64,Lx::Int64,Ly::Int64,Ns::StepRange{Int64, Int64})
    
    for N in Ns

        ops = ham(U,Lx,Ly) # create the Hamiltonian
        nup = N%2 == 0 ? N ÷ 2 : (N+1) ÷ 2
        ndn = N%2 == 0 ? N ÷ 2 : (N-1) ÷ 2

        block = Electron(Lx*Ly, nup, ndn)   # create Hilbert space
        res = eigs_lanczos(ops, block, neigvals = 3)
        eigs = res.eigenvalues

        filename = "data/tos_ahm/U($U)_N($N)_Lx($Lx)_Ly($Ly).h5"

        h5open(filename,"w") do f
            write(f,"spectrum",eigs)
        end

    end

end

Us = [0.0,-2.0,-10.0]
Lx, Ly = (4,4)
Nmin = 4
Nmax = 2*Lx*Ly - 4
Ns = Nmin:2:Nmax

say_hello()
for U in Us
    println("Now doing U=$U")
    main(U,Lx,Ly,Ns)
end
