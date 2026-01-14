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

#--8<-- [start:usage_guide_spm1]
coo_mat = coo_matrix(ops, block)
csr_mat = csr_matrix(ops, block)
csc_mat = csc_matrix(ops, block)
# --8<-- [end:usage_guide_spm1]

#--8<-- [start:usage_guide_spm2]
using SparseArrays
csc_mat = csc_matrix(ops, block)
jmat = SparseMatrixCSC(csc_mat.nrows, csc_mat.ncols, csc_mat.colptr, csc_mat.row, csc_mat.data)
# --8<-- [end:usage_guide_spm2]

#--8<-- [start:usage_guide_spm3]
# Heisenberg chain with open boundary conditions
N = 22
block = Spinhalf(N)
ops = OpSum()
for i in 1:(N-1)
    ops += "J" * Op("SdotS", [i, i+1])
end
ops["J"] = 1.0

# perform default (matrix-free) Lanczos iterations
@time begin
    println(eigval0(ops, block))
end

# perform Lanczos iterations on sparse CSR Hamiltonian
@time begin
    csr_mat = csr_matrix(ops, block)
    println(eigval0(csr_mat, block))
end

# Example output:
#
#   -9.568075875976074
#       6.859141 seconds (7 allocations: 512 bytes)
#   -9.568075875976076
#       5.954734 seconds (376 allocations: 2.000 GiB, 16.86% gc time)
#--8<-- [end:usage_guide_spm3]

    
end
