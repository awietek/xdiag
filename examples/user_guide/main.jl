using XDiag
using LinearAlgebra

let
# --8<-- [start:usage_guide_hs1]
N = 8
hs = Spinhalf(8)
# --8<-- [end:usage_guide_hs1]

# --8<-- [start:usage_guide_hs2]
for spins in hs
    println(to_string(spins))
    println(index(hs, spins))
end
println("dim: ", size(hs))
# --8<-- [end:usage_guide_hs2]

# --8<-- [start:usage_guide_hs3]
nup = 2
b1 = Spinhalf(N, nup)

ndn = 1
b2 = tJ(N, nup, ndn)
b3 = Electron(N, nup, ndn)
# --8<-- [end:usage_guide_hs3]

# --8<-- [start:usage_guide_op1]
ops = OpSum()
for i in 1:N
    s1 = i
    s2 = mod1(i+1, N)
    ops += "J" * Op("SdotS", [s1, s2])
end
ops["J"] = 1.0
# --8<-- [end:usage_guide_op1]

# --8<-- [start:usage_guide_mat1]
H = matrix(ops, block);
display(H)
# --8<-- [end:usage_guide_mat1]

# --8<-- [start:usage_guide_mat2]
(evals, evecs) = eigen(Symmetric(H))
# --8<-- [end:usage_guide_mat2]

# --8<-- [start:usage_guide_stat1]
real = true
psi1 = State(b, real=real)
psi2 = zero_state(b, real=real)
# --8<-- [end:usage_guide_stat1]

# --8<-- [start:usage_guide_stat2]
d = size(block)
v = rand(d)
psi = State(block, v)
# --8<-- [end:usage_guide_stat2]

# --8<-- [start:usage_guide_stat3]
psi1 = product_state(block, ["Up", "Dn"])
psi2 = random_state(block)
# --8<-- [end:usage_guide_stat3]

# --8<-- [start:usage_guide_stat4]
nrm = norm(psi)
dot = dot(psi1, psi2)
# --8<-- [end:usage_guide_stat4]

# --8<-- [start:usage_guide_stat5]
v = vector(psi)
# --8<-- [end:usage_guide_stat5]

# --8<-- [start:usage_guide_stat6]
phi = apply(H, psi)
# --8<-- [end:usage_guide_stat6]

# --8<-- [start:usage_guide_iter1]
e0 = eigval0(H, block)
# --8<-- [end:usage_guide_iter1]

# --8<-- [start:usage_guide_iter2]
e0, psi0 = eig0(H, block)
# --8<-- [end:usage_guide_iter2]

# --8<-- [start:usage_guide_iter3]
t = 1.0
phi = time_evolve(H, psi0, t)
# --8<-- [end:usage_guide_iter3]

# --8<-- [start:usage_guide_measu1]
for i in 1:N
    op = Op("SzSz", {1, i})
    corr = inner(op, psi0)
end
# --8<-- [end:usage_guide_measu1]

# --8<-- [start:usage_guide_io1]
fl = FileToml("spinhalf_chain.toml")
ops = read_opsum(fl, "Interactions")
ops["J"] = 1.0
# --8<-- [end:usage_guide_io1]

# --8<-- [start:usage_guide_sym1]
T = Permutation([2, 3, 4, 5, 6, 7, 8, 1])
# --8<-- [end:usage_guide_sym1]

# --8<-- [start:usage_guide_sym2]
group = PermutationGroup([T^0, T^1, T^2, T^3, T^4, T^5, T^6, T^7])
# --8<-- [end:usage_guide_sym2]

# --8<-- [start:usage_guide_sym3]
chi = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
k = Representation(group, chi)
# --8<-- [end:usage_guide_sym3]

# --8<-- [start:usage_guide_sym4]
blk = Spinhalf(N, nup, irrep)
# --8<-- [end:usage_guide_sym4]

# --8<-- [start:usage_guide_sym5]
for spins in blk
    println(to_string(spins))
end
# --8<-- [end:usage_guide_sym5]

# --8<-- [start:usage_guide_sym6]
fl = FileToml("symmetries.toml")
group = read_permutation_group(fl, "Symmetries")
irrep = read_representation(fl, "k.zero", "Symmetries")
# --8<-- [end:usage_guide_sym6]

#--8<-- [start:usage_guide_sym7]
og = symmetrize(ops, group)
oi = symmetrize(ops, irrep)
# --8<-- [end:usage_guide_sym7]












# --8<-- [start:usage_guide_5]
ops = OpSum()
for i in 1:N
    ops += "J" * Op("SdotS", [i, mod1(i+1, N)])
end
ops["J"] = 1.0
# --8<-- [end:usage_guide_5]

# --8<-- [start:usage_guide_6]
e0, psi0 = eig0(ops, block)
println("e0:", e0)
# --8<-- [end:usage_guide_6]

# --8<-- [start:usage_guide_7]
H = matrix(ops, block); 
(evals, evecs) = eigen(Symmetric(H))
println("e0:", evals[1], "e1:", evals[2]);
# --8<-- [end:usage_guide_7]

# --8<-- [start:usage_guide_8]
for i in 2:N
    op = Op("SzSz", [1, i])
    corr = inner(op, psi0)
    println("<Sz_0 Sz_$i> = ", corr)
end
# --8<-- [end:usage_guide_8]

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
