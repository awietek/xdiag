using XDiag
using Printf
using HDF5

function SpectrumSz()
    energies = Vector{Float64}[] ## Collect the energies
    IrrepList = Vector{String}[] ## Collect the energies
    Nsites = 18
    numbeig = 6 # number of lanczos vectors to converge

    Irrep = []

    fl = FileToml("triangular.18.10418.J1J2.sublattices.tsl.toml") #TOML file with Shastry-Sutherland Interactions 
    ops = read_opsum(fl, "Interactions")
    # Define couplings
    ops["J1"] = 1.0
    ops["J2"] = 0.0 #  For this ratio of couplings we are in the 120 phase.

    Irreps =  ["Gamma.C2.A", "Gamma.C2.B", "K0.C1.A", "K1.C1.A", "M.C2.A", "M.C2.B", "X0.C1.A","X1.C1.A","X2.C1.A", "Z0.C1.A","Z1.C1.A","0.C1.A","1.C1.A","2.C1.A"]

    # Loop different Irreps
    for j in Irreps
        irrep = read_representation(fl, j)
        # Loop different total magnetization sector:
        for nup in 0:Nsites
            block = Spinhalf(Nsites, nup, irrep)
            r = eigvals_lanczos(ops, block, neigvals=numbeig)
            eig = r.eigenvalues
            for i in 1:length(eig)
                push!(energies, [nup, eig[i]])
                push!(IrrepList,[j])
            end
        end
    end

    filename = @sprintf("energies_tower_of_states.triangular.Nsites.%d.outfile.h5", Nsites)
    h5open(filename, "w") do file
        write(file, "energies", hcat(energies...))
        write(file, "irrep", hcat(IrrepList...))
    end
end

SpectrumSz()