using XDiag

# --8<-- [start:Permutation]
p1 = Permutation([0, 2, 1, 3])
p2 = Permutation([2, 0, 1, 3])

@show inverse(p1)
@show p1 * p2
# --8<-- [end:Permutation]

# --8<-- [start:PermutationGroup]
# Define a cyclic group of order 3
p1 = Permutation([0, 1, 2])
p2 = Permutation([1, 2, 0])
p3 = Permutation([2, 0, 1])
C3 = PermutationGroup([p1, p2, p3])

@show size(C3)
@show n_sites(C3)
@show inverse(C3, 1) # = 2
# --8<-- [end:PermutationGroup]


# --8<-- [start:Representation]
r1 = Representation([1, -1, 1, -1])
r2 = Representation([1, 1im, -1, -1im])

@show r1 * r2
# --8<-- [end:Representation]

# --8<-- [start:Spinhalf]
N = 4
nup = 2

# without Sz conservation
block = Spinhalf(N)
@show block


# with Sz conservation
block_sz = Spinhalf(N, nup)
@show block_sz

# with symmetries, without Sz
p1 = Permutation([0, 1, 2, 3])
p2 = Permutation([1, 2, 3, 0])
p3 = Permutation([2, 3, 0, 1])
p4 = Permutation([3, 0, 1, 2])
group = PermutationGroup([p1, p2, p3, p4])
rep = Representation([1, -1, 1, -1])
block_sym = Spinhalf(N, group, rep)
@show block_sym

# with symmetries and Sz
block_sym_sz = Spinhalf(N, nup, group, rep)
@show block_sym_sz

@show n_sites(block_sym_sz)
@show size(block_sym_sz)

# Iteration
for pstate in block_sym_sz
    @show pstate, index(block_sym_sz, pstate)
end
@show permutation_group(block_sym_sz)
@show irrep(block_sym_sz)
# --8<-- [end:Spinhalf]

# --8<-- [start:tJ]
N = 4
nup = 2
ndn = 1

# without permutation symmetries
block = tJ(N, nup, ndn)
@show block

# with permutation symmetries
p1 = Permutation([0, 1, 2, 3])
p2 = Permutation([1, 2, 3, 0])
p3 = Permutation([2, 3, 0, 1])
p4 = Permutation([3, 0, 1, 2])
group = PermutationGroup([p1, p2, p3, p4])
rep = Representation([1, -1, 1, -1])
block_sym = tJ(N, nup, ndn, group, rep)
@show block_sym

@show n_sites(block_sym)
@show size(block_sym)

# Iteration
for pstate in block_sym
    @show pstate, index(block_sym, pstate)
end
@show permutation_group(block_sym)
@show irrep(block_sym)
# --8<-- [end:tJ]


# --8<-- [start:Electron]
N = 4
nup = 2
ndn = 1

# without number conservation
block = Electron(N)
@show block

# with number conservation
block_np = Electron(N, nup, ndn)
@show block_np

# with symmetries, without number conservation
p1 = Permutation([0, 1, 2, 3])
p2 = Permutation([1, 2, 3, 0])
p3 = Permutation([2, 3, 0, 1])
p4 = Permutation([3, 0, 1, 2])
group = PermutationGroup([p1, p2, p3, p4])
rep = Representation([1, -1, 1, -1])
block_sym = Electron(N, group, rep)
@show block_sym

# with symmetries and number conservation
block_sym_np = Electron(N, nup, ndn, group, rep)
@show block_sym_np

@show n_sites(block_sym_np)
@show size(block_sym_np)

# Iteration
for pstate in block_sym_np
    @show pstate, index(block_sym_np, pstate)
end
@show permutation_group(block_sym_np)
@show irrep(block_sym_np)
# --8<-- [end:Electron]


# --8<-- [start:matrix]
let
    # Creates matrix H_{k=2} in Eq (18.23) of https://link.springer.com/content/pdf/10.1007/978-3-540-74686-7_18.pdf
    N = 4
    nup = 3
    ndn = 2

    # Define a Hubbard chain model
    ops = OpSum()
    for i in 1:N
        ops += Op("HOP", "T", [i-1, i % N])
    end
    ops["T"] = 1.0;
    ops["U"] = 5.0;

    # Create the a permutation group
    p1 = Permutation([0, 1, 2, 3])
    p2 = Permutation([1, 2, 3, 0])
    p3 = Permutation([2, 3, 0, 1])
    p4 = Permutation([3, 0, 1, 2])
    group = PermutationGroup([p1, p2, p3, p4])
    irrep = Representation([1, -1, 1, -1])
    block = Electron(N, nup, ndn, group, irrep)

    H = matrix(ops, block)
    display(H)
end
# --8<-- [end:matrix]

# --8<-- [start:eigval0]
let 
    N = 8;
    nup = N ÷ 2;
    block = Spinhalf(N, nup);
    
    # Define the nearest-neighbor Heisenberg model
    ops = OpSum()
    for i in 1:N
        ops += Op("HB", "J", [i-1, i % N])
    end
    ops["J"] = 1.0;

    e0 = eigval0(ops, block);
end
# --8<-- [end:eigval0]

# --8<-- [start:eig0]
let 
    N = 8;
    nup = N ÷ 2;
    block = Spinhalf(N, nup);
    
    # Define the nearest-neighbor Heisenberg model
    ops = OpSum()
    for i in 1:N
        ops += Op("HB", "J", [i-1, i % N])
    end
    ops["J"] = 1.0;

    e0, gs = eig0(ops, block);
end
# --8<-- [end:eig0]


# --8<-- [start:op]
op = Op("HOP", "T", [1, 2])
@show op
@show type(op)
@show convert(String, coupling(op))
@show size(op), op[1], op[2]
@show sites(op) == [1, 2]
@show isexplicit(op)

op = Op("HOP", 1.23, [1, 2])
@show op
@show isreal(op)
@show ismatrix(op)
@show isexplicit(op)

op = Op("SY", [0 -im; im 0], 1)
@show op
@show isreal(op)
@show ismatrix(op)
@show isexplicit(op)
# --8<-- [end:op]


# --8<-- [start:opsum]
# Define the 1D transverse-field Ising chain
let 
    N = 12
    J = 1.0
    h = 0.5
    Sx = [0 1; 1 0]

    ops = OpSum()
    for i in 1:N
        ops += Op("ISING", "J", [i-1, i % N])
        ops += Op("SX", h * Sx, i-1)
    end

    ops["J"] = 1.0;
    @show ops
    @show defined(ops, "J")
    @show isreal(ops)
    @show isexplicit(ops)
end
# --8<-- [end:opsum]

# --8<-- [start:coupling]
cpl = Coupling("J")
@show type(cpl)
@show isexplicit(cpl)

cpl = Coupling(1.23)
@show ismatrix(cpl)
@show convert(Float64, cpl)
@show convert(ComplexF64, cpl)

cpl = Coupling([1 2; -2 1])
@show ismatrix(cpl)
@show isreal(cpl)
@show convert(Matrix{Float64}, cpl)
@show convert(Matrix{ComplexF64}, cpl)
# --8<-- [end:coupling]

# --8<-- [start:state]
block = Spinhalf(2)
psi1 = State(block, [1.0, 2.0, 3.0, 4.0])
@show psi1
display(vector(psi1))
make_complex!(psi1)
display(vector(psi1))

psi2 = State(block, real=false, n_cols=3)
@show psi2
display(matrix(psi2))

psi3 = State(block, [1.0+4.0im, 2.0+3.0im, 3.0+2.0im, 4.0+1.0im])
display(vector(psi3))
display(vector(real(psi3)))
display(vector(imag(psi3)))
# --8<-- [end:state]


# --8<-- [start:product_state]
pstate = ProductState(["Up", "Dn", "Emp", "UpDn"])
for s in pstate
    @show s
end
@show pstate

pstate = ProductState()
push!(pstate, "Dn")
push!(pstate, "Up")
push!(pstate, "Dn")
@show n_sites(pstate)
for s in pstate
    @show s
end
@show pstate
# --8<-- [end:product_state]
