using HDF5, CairoMakie

# === Load Data ===
# Open the HDF5 file and read datasets "energies" and "irrep"
data18 = h5open("energies_tower_of_states.triangular.Nsites.18.outfile.h5", "r")
AllEnergies = read(data18["energies"])
Allirreps = read(data18["irrep"])[1,:]
close(data18)

# === Function Definitions ===

# Computes total spin values given energies and irreps.
function get_STotal(energies, irreps, Nspins)
    # In Julia the columns are 1-indexed; here column 2 corresponds to Python’s energies[:,1]
    energies2 = round.(energies[2, :], digits=8)
    E0 = unique(energies2)
    sort!(E0)
    Stot = similar(E0)
    Irreps_arr = Vector{Any}(undef, length(E0))

    for (i, e0) in enumerate(E0)
        mask = findall(x -> x == e0, energies2)
        # Column 1 of energies corresponds to Python’s energies[:,0]
        Sz = energies[1, mask] .- Nspins / 2
        vals = Sz .* (Sz .+ 1)
        max_val, arg = findmax(vals)
        Stot[i] =abs(max_val)
        # Get the corresponding irreps element.
        Irreps_arr[i] = irreps[mask][arg]
    end
    return Stot, E0, Irreps_arr
end

# Extracts lower energy levels for a specified irreducible representation.
function get_lower(Stotal, sortedEnergies, allenergies, IrrepList, IrreListAll, irr)
    target = irr[1]
    maskAll = findall(x -> x == target, IrreListAll)

    mask = findall(x -> x == target, IrrepList)

    energies2 = round.(allenergies[2, maskAll], digits=8)
    

    St = round.(Stotal, digits=0)
    # Consider only the entries corresponding to the target irreducible rep.
    St_mask = St[mask]
    
    St2 = sort(unique(St_mask))
    
    E2 = sortedEnergies[mask]
    
    Ef = zeros(length(St2), 2)

    for (i, s) in enumerate(St2)
        mask2 = findall(x -> x == s, St_mask)

        E_subset = E2[mask2]
        
        Ef[i, 1] = minimum(E_subset)
        
        multiplicity = findall(x -> x == Ef[i, 1], energies2)
        
        Ef[i, 2] = length(multiplicity)
    end
    return St2, Ef
end



function make_Plot()
    # --- Compute Quantities ---
    Nspins = 18
    Stot, SEtot, SortedIr = get_STotal(AllEnergies, Allirreps,Nspins)

    SGamma, EGamma = get_lower(Stot, SEtot, AllEnergies, SortedIr, Allirreps, ["Gamma.C2.A"])
    SGammaB, EGammaB = get_lower(Stot, SEtot, AllEnergies, SortedIr, Allirreps, ["Gamma.C2.B"])
    SGammaK, EGammaK = get_lower(Stot, SEtot, AllEnergies, SortedIr, Allirreps, ["K0.C1.A"])

    # === Plotting ===
    fig = Figure(ratio=1.618)
    ax = Axis(fig[1, 1],
        xlabel=L"S_{\text{total}}(S_{\text{ total}}+1)",
        ylabel=L"$(E - E_{\text{GS} })/J_1$" )

    xlims!(ax, -0.5, 12.5),
    ylims!(ax, -0.4, 5)
    # Plot the reference line: Stot vs. 0.265*Stot
    lines!(ax, Stot, 0.265 .* Stot, color=:black)

    Egs = minimum(SEtot)
    # Plot all energies shifted by the ground state energy.
    scatter!(ax, Stot, SEtot .- Egs, color=:black, markersize=4)

    maskA = [2, 4]  

    print(SGamma)

    scatter!(ax, SGamma[maskA], EGamma[maskA, 1] .- Egs,
        marker=:utriangle, markersize=10,
        label=L"$\Gamma.\mathrm{C2.A}$")
    @info "Multiplicity for Gamma.C2.A: $(EGamma[maskA, 2])"

    maskB = [1, 3, 4]  
    scatter!(ax, SGammaB[maskB], EGammaB[maskB, 1] .- Egs,
        marker=:utriangle, markersize=10,
        label=L"$\Gamma.\mathrm{C2.B}$")
    @info "Multiplicity for Gamma.C2.B: $(EGammaB[maskB, 2])"

    maskC = [2, 3, 4]  
    scatter!(ax, SGammaK[maskC], EGammaK[maskC, 1] .- Egs,
        marker=:star5, markersize=10,
        label=L"$K_0.\mathrm{C1.A}$")

    @info "Multiplicity for K0.C1.A: $(EGammaK[maskC, 2])"

    text!(ax, L"J_2=0", position=(10, 4.5), align=(:left, :center))


    axislegend(ax; framevisible=true, position=:rb, labelsize=15)
    display(fig)

    return nothing
end


make_Plot()
