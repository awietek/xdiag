using XDiag
using LinearAlgebra

let
# --8<-- [start:first_steps_1]
N = 8
hspace = Spinhalf(8)
# --8<-- [end:first_steps_1]

# --8<-- [start:first_steps_2]
for spins in hspace
    println(to_string(spins))
end
# --8<-- [end:first_steps_2]

# --8<-- [start:first_steps_3]
nup = 4;
block = Spinhalf(N, nup);
for spins in block
    println(to_string(spins))
end
# --8<-- [end:first_steps_3]

# --8<-- [start:first_steps_4]
@show size(hspace)
@show size(block)
# --8<-- [end:first_steps_4]

# --8<-- [start:first_steps_5]
ops = OpSum()
for i in 1:N
    ops += "J" * Op("SdotS", [i, mod1(i+1, N)])
end
ops["J"] = 1.0
# --8<-- [end:first_steps_5]

# --8<-- [start:first_steps_6]
e0, psi0 = eig0(ops, block)
println("e0:", e0)
# --8<-- [end:first_steps_6]

# --8<-- [start:first_steps_7]
H = matrix(ops, block); 
(evals, evecs) = eigen(Symmetric(H))
println("e0:", evals[1], "e1:", evals[2]);
# --8<-- [end:first_steps_7]

# --8<-- [start:first_steps_8]
for i in 2:N
    op = Op("SzSz", [1, i])
    corr = inner(op, psi0)
    println("<Sz_0 Sz_$i> = ", corr)
end
# --8<-- [end:first_steps_8]

# --8<-- [start:io_1]
# fl = FileToml("/Users/awietek/Research/Software/xdiag/examples/user_guide/spinhalf_chain.toml")
# ops = read_opsum(fl, "Interactions")
# --8<-- [end:io_1]


# --8<-- [start:symmetries_1]
T = Permutation([2, 3, 4, 5, 6, 7, 8, 1])
# --8<-- [end:symmetries_1]
 

# --8<-- [start:symmetries_2]
group = PermutationGroup([T^0, T^1, T^2, T^3, T^4, T^5, T^6, T^7])
# --8<-- [end:symmetries_2]

# --8<-- [start:symmetries_3]
irrep_k_0 = Representation(group, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
irrep_k_pi = Representation(group, [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0])
# --8<-- [end:symmetries_3]

# --8<-- [start:symmetries_4]
block_k_0 = Spinhalf(N, nup, irrep_k_0)
block_k_pi = Spinhalf(N, nup, irrep_k_pi)
e0_k_0 = eigval0(ops, block_k_0)
e0_k_pi = eigval0(ops, block_k_pi)
println("e0: k=0: $e0_k_0, k=pi: $e0_k_pi")
# --8<-- [end:symmetries_4]

    
end
