using XDiag

# --8<-- [start:Permutation]
p1 = Permutation([1, 3, 2, 4])
p2 = Permutation([3, 1, 2, 4])

@show inv(p1)
@show p1 * p2
@show p1 ^ 2
# --8<-- [end:Permutation]

# --8<-- [start:PermutationGroup]
# Define a cyclic group of order 3
p1 = Permutation([1, 2, 3])
p2 = Permutation([2, 3, 1])
p3 = Permutation([3, 1, 2])
C3 = PermutationGroup([p1, p2, p3])

@show size(C3)
@show nsites(C3)
# --8<-- [end:PermutationGroup]


# --8<-- [start:Representation]
p = Permutation([2, 3, 4, 1])
C4 = PermutationGroup([p^0, p^1, p^2, p^3])
r1 = Representation(C4, [1.0, -1.0, 1.0, -1.0])
r2 = Representation(C4, [1.0, 1.0im, -1.0, -1.0im])
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
p = Permutation([2, 3, 4, 1])
group = PermutationGroup([p^0, p^1, p^2, p^3])
rep = Representation(group, [1.0, -1.0, 1.0, -1.0])
block_sym = Spinhalf(N, rep)
@show block_sym

# with symmetries and Sz
block_sym_sz = Spinhalf(N, nup, rep)
@show block_sym_sz
@show nsites(block_sym_sz)
@show size(block_sym_sz)

# Iteration
for pstate in block_sym_sz
    @show pstate, index(block_sym_sz, pstate)
end
# --8<-- [end:Spinhalf]

# --8<-- [start:tJ]
N = 4
nup = 2
ndn = 1

# without permutation symmetries
block = tJ(N, nup, ndn)
@show block

# with permutation symmetries
p = Permutation([2, 3, 4, 1])
group = PermutationGroup([p^0, p^1, p^2, p^3])
rep = Representation(group, [1.0, -1.0, 1.0, -1.0])
block_sym = tJ(N, nup, ndn, rep)
@show block_sym

@show nsites(block_sym)
@show size(block_sym)

# Iteration
for pstate in block_sym
    @show pstate, index(block_sym, pstate)
end
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
p = Permutation([2, 3, 4, 1])
group = PermutationGroup([p^0, p^1, p^2, p^3])
rep = Representation(group, [1.0, -1.0, 1.0, -1.0])
block_sym = Electron(N, rep)
@show block_sym

# with symmetries and number conservation
block_sym_np = Electron(N, nup, ndn, rep)
@show block_sym_np
@show nsites(block_sym_np)
@show size(block_sym_np)

# Iteration
for pstate in block_sym_np
    @show pstate, index(block_sym_np, pstate)
end
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
        ops += "T" * Op("Hop", [i, mod1(i+1, N)])
    end
    ops += "U" * Op("HubbardU")
    ops["T"] = 1.0;
    ops["U"] = 5.0;

    # Create the a permutation group
    p = Permutation([2, 3, 4, 1])
    group = PermutationGroup([p^0, p^1, p^2, p^3])
    irrep = Representation(group, [1.0, -1.0, 1.0, -1.0])
    block = Electron(N, nup, ndn, irrep)

    H = matrix(ops, block)
    display(H)
end
# --8<-- [end:matrix]


# --8<-- [start:coo_matrix]
let 
    N = 8
    nup = N ÷ 2
    block = Spinhalf(N, nup)
    
    # Define the nearest-neighbor Heisenberg model
    ops = OpSum()
    for i in 1:N
        ops += "J" * Op("SdotS", [i, mod1(i+1, N)])
    end
    ops["J"] = 1.0
    spmat = coo_matrix(ops, block)
    spmat_32 = coo_matrix_32(ops, block)
end
# --8<-- [end:coo_matrix]


# --8<-- [start:csc_matrix]
let 
    N = 8
    nup = N ÷ 2
    block = Spinhalf(N, nup)
    
    # Define the nearest-neighbor Heisenberg model
    ops = OpSum()
    for i in 1:N
        ops += "J" * Op("SdotS", [i, mod1(i+1, N)])
    end
    ops["J"] = 1.0
    spmat = csc_matrix(ops, block)
    spmat_32 = csc_matrix_32(ops, block)
end
# --8<-- [end:csc_matrix]

# --8<-- [start:csr_matrix]
let 
    N = 8
    nup = N ÷ 2
    block = Spinhalf(N, nup)
    
    # Define the nearest-neighbor Heisenberg model
    ops = OpSum()
    for i in 1:N
        ops += "J" * Op("SdotS", [i, mod1(i+1, N)])
    end
    ops["J"] = 1.0
    spmat = csr_matrix(ops, block)
    spmat_32 = csr_matrix_32(ops, block)
end
# --8<-- [end:csr_matrix]

# --8<-- [start:eigval0]
let 
    N = 8
    nup = N ÷ 2
    block = Spinhalf(N, nup)
    
    # Define the nearest-neighbor Heisenberg model
    ops = OpSum()
    for i in 1:N
        ops += "J" * Op("SdotS", [i, mod1(i+1, N)])
    end
    ops["J"] = 1.0
    e0 = eigval0(ops, block);
end
# --8<-- [end:eigval0]

# --8<-- [start:eig0]
let 
    N = 8
    nup = N ÷ 2
    block = Spinhalf(N, nup)
    
    # Define the nearest-neighbor Heisenberg model
    ops = OpSum()
    for i in 1:N
        ops += "J" * Op("SdotS", [i, mod1(i+1, N)])
    end
    ops["J"] = 1.0
    e0, gs = eig0(ops, block);
end
# --8<-- [end:eig0]

# --8<-- [start:eigvals_lanczos]
let
    N = 8
    nup = N ÷ 2
    block = Spinhalf(N, nup)
    
    # Define the nearest-neighbor Heisenberg model
    ops = OpSum()
    for i in 1:N
        ops += "J" * Op("SdotS", [i, mod1(i+1, N)])
    end
    ops["J"] = 1.0

    # With random intial state
    res = eigvals_lanczos(ops, block)
    @show res.alphas
    @show res.betas
    @show res.eigenvalues

    # With specific initial state
    psi0 = product_state(block, ["Up", "Dn", "Up", "Dn", "Up", "Dn", "Up", "Dn"])
    res2 = eigvals_lanczos(ops, psi0)
    @show res.alphas
    @show res.betas
    @show res.eigenvalues
end
# --8<-- [end:eigvals_lanczos]


# --8<-- [start:eigs_lanczos]
let
    N = 8
    nup = N ÷ 2
    block = Spinhalf(N, nup)
    
    # Define the nearest-neighbor Heisenberg model
    ops = OpSum()
    for i in 1:N
        ops += "J" * Op("SdotS", [i, mod1(i+1, N)])
    end
    ops["J"] = 1.0

    # With random intial state
    res = eigs_lanczos(ops, block)
    @show res.alphas
    @show res.betas
    @show res.eigenvalues
    @show res.eigenvectors
end
# --8<-- [end:eigs_lanczos]


# --8<-- [start:time_evolve]
let
    N = 8
    nup = N ÷ 2
    block = Spinhalf(N, nup)
    
    # Define the nearest-neighbor Heisenberg model
    ops = OpSum()
    for i in 1:N
        ops += "J" * Op("SdotS", [i, mod1(i+1, N)])
    end
    ops["J"] = 1.0

    psi0 = product_state(block, ["Up", "Dn", "Up", "Dn", "Up", "Dn", "Up", "Dn"])
    time = 1.0
    psi = time_evolve(ops, psi0, time)
    time_evolve_inplace(ops, psi0, time)
    @show isapprox(psi0, psi)
end
# --8<-- [end:time_evolve]

# --8<-- [start:imaginary_time_evolve]
let
    N = 8
    nup = N ÷ 2
    block = Spinhalf(N, nup)
   
    # Define the nearest-neighbor Heisenberg model
    ops = OpSum()
    for i in 1:N
        ops += Op("SdotS", [i, mod1(i+1, N)])
    end

    # Compute ground state energy
    e0 = eigval0(ops, block)
 
    psi0 = product_state(block, ["Up", "Dn", "Up", "Dn", "Up", "Dn", "Up", "Dn"])
    time = 1.0
    psi = imaginary_time_evolve(ops, psi0, time,
                                precision=1e-12, shift=e0)
    imaginary_time_evolve_inplace(ops, psi0, time, precision=1e-12, shift=e0)
    @show isapprox(psi0, psi)
end
# --8<-- [end:imaginary_time_evolve]

# --8<-- [start:evolve_lanczos]
let
    N = 8
    nup = N ÷ 2
    block = Spinhalf(N, nup)
    
    # Define the nearest-neighbor Heisenberg model
    ops = OpSum()
    for i in 1:N
        ops += Op("SdotS", [i, mod1(i+1, N)])
    end

    # Compute ground state energy
    e0 = eigval0(ops, block)
 
    psi0 = product_state(block, ["Up", "Dn", "Up", "Dn", "Up", "Dn", "Up", "Dn"])
    time = 1.0
    res = evolve_lanczos(ops, psi0, time, precision=1e-12, shift=e0, normalize=true)
    @show res.alphas
    @show res.betas
end
# --8<-- [end:evolve_lanczos]

    
# --8<-- [start:time_evolve_expokit]
let
    N = 8
    nup = N ÷ 2
    block = Spinhalf(N, nup)
    
    # Define the nearest-neighbor Heisenberg model
    ops = OpSum()
    for i in 1:N
        ops += Op("SdotS", [i, mod1(i+1, N)])
    end

    psi0 = product_state(block, ["Up", "Dn", "Up", "Dn", "Up", "Dn", "Up", "Dn"])
    time = 1.0
    res1 = time_evolve_expokit(ops, psi0, time, precision=1e-8)
    res2 = time_evolve_expokit_inplace(ops, psi0, time, precision=1e-8)
    @show isapprox(psi0, res1.state)
    @show res1.error
    @show res1.hump
end
# --8<-- [end:time_evolve_expokit]


# --8<-- [start:Op]
op = "T" * Op("Hop", [1, 2])
@show op

op = 1.23 * Op("Hop",  [1, 2])
@show op

op = Op("Matrix", 1, [0 -1.0im; 1.0im 0])
@show op
@show isreal(op)
# --8<-- [end:Op]


# --8<-- [start:opsum]
# Define the 1D transverse-field Ising chain
let 
    N = 12
    J = 1.0
    h = 0.5
    Sx = [0.0 1.0; 1.0 0.0]

    ops1 = OpSum()
    for i in 1:N
        ops1 += J * Op("SzSz", [i, mod1(i+1, N)])
        ops1 += h * Op("Matrix", i, Sx)
    end

    ops2 = OpSum()
    for i in 1:N
        ops2 += "J" * Op("SzSz", [i, mod1(i+1, N)])
        ops2+= "h" * Op("Matrix", i, Sx)
    end
    ops2["J"] = J;
    ops2["h"] = h;
    
    @show isapprox(ops1, ops2)
    @show isapprox(ops1 + ops2, 2.0 * ops1)
    @show to_string(ops1)
end
# --8<-- [end:opsum]


# --8<-- [start:hc]
cdagup = Op("Cdagup", 1)
sdots = Op("SdotS", [1, 2])
hop = (1.0 + 1.0im) * Op("Hop", [1, 2])
@show cdagup == hc(cdagup)
@show sdots == hc(sdots)
@show hop == hc(hop)
# --8<-- [end:hc]


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
@show nsites(pstate)
for s in pstate
    @show s
end
@show pstate
# --8<-- [end:product_state]


# --8<-- [start:random_state]
block = Spinhalf(2)
state = State(block, real=false)  # complex State
rstate1 = RandomState(1234)
fill(state, rstate1)
display(vector(state))

rstate2 = RandomState(4321)
fill(state, rstate2)
display(vector(state))

fill(state, rstate1)
display(vector(state))
# --8<-- [end:random_state]

# --8<-- [start:fill]
block = Spinhalf(2)
state = State(block)
pstate = ProductState(["Up", "Dn"])
fill(state, pstate)
display(vector(state))

rstate = RandomState(1234)
fill(state, rstate)
display(vector(state))
# --8<-- [end:fill]


# --8<-- [start:create_state]
block = Spinhalf(2)
state = product_state(block, ["Up", "Dn"])
display(vector(state))

zero(state)
display(vector(state))

state = random_state(block, real=false, seed=1234, normalized=true)
display(vector(state))

state = zero_state(block, real=true, ncols=2)
display(matrix(state))
# --8<-- [end:create_state]


# --8<-- [start:algebra]
let 
    N = 8
    block = Spinhalf(N,  N ÷ 2)
    ops = OpSum()
    for i in 1:N
        ops += Op("SdotS", [i, mod1(i+1, N)])
    end
    e0, psi = eig0(ops, block);

    @show norm(psi)
    @show norm1(psi)
    @show norminf(psi)
    @show dot(psi, psi)
    @show e0, inner(ops, psi)

    phi = random_state(block)
    display(vector(phi))
    display(vector(psi))
    display(vector(psi + 2.0*phi))
    display(vector(psi*3.0im + phi/2.0))
end
# --8<-- [end:algebra]

# --8<-- [start:apply]
let 
    N = 8
    block = Spinhalf(N,  N ÷ 2)
    ops = OpSum()
    for i in 1:N
        ops += Op("SdotS", [i, mod1(i+1, N)])
    end
    e0, psi = eig0(ops, block);
    phi = apply(Op("S+", 2), psi)
    @show inner(ops, psi)
    @show inner(ops, phi)
end
# --8<-- [end:apply]


# --8<-- [start:sparse_apply]
let
    # Real vector apply
    N = 8
    block = Spinhalf(N,  N ÷ 2)
    ops = OpSum()
    for i in 1:N
        ops += Op("SdotS", [i, mod1(i+1, N)])
    end
    spmat = csr_matrix(ops, block)
    rstate = random_state(block)
    res = apply(spmat, vector(rstate))

    # Complex matrix apply
    N = 4
    block = Electron(N,  N ÷ 2, N ÷ 2)
    ops = OpSum()
    for i in 1:N
        ops += (1.0 + 1.0im) * Op("Hop", [i, mod1(i+1, N)])
    end
    spmat = csr_matrix(ops, block)
    ncols = 3
    rstate = random_state(block; real=false, ncols=ncols)
    in = matrix(rstate)
    out = zeros(ComplexF64, size(block), ncols)
    apply(spmat, in, out)
end
# --8<-- [end:sparse_apply]

# --8<-- [start:symmetrize]
let
    N = 4
    nup = 2
    block = Spinhalf(N, nup)
    p1 = Permutation([1, 2, 3, 4])
    p2 = Permutation([2, 3, 4, 1])
    p3 = Permutation([3, 4, 1, 2])
    p4 = Permutation([4, 1, 2, 3])
    group = PermutationGroup([p1, p2, p3, p4])
    rep = Representation(group)
    block_sym = Spinhalf(N, rep)

    ops = OpSum()
    for i in 1:N
        ops += Op("SdotS", [i, mod1(i+1, N)])
    end

    e0, psi = eig0(ops, block);
    e0, psi_sym = eig0(ops, block_sym);

    corr = Op("SdotS", [1, 2])
    nn_corr = inner(corr, psi)
    corr_sym = symmetrize(corr, group)
    nn_corr_sym = inner(corr_sym, psi_sym)
    @show nn_corr, nn_corr_sym
end
# --8<-- [end:symmetrize]


# --8<-- [start:read_opsum]
file = "triangular.9.hop.sublattices.tsl.toml"
fl = FileToml(file)
ops = read_opsum(fl, "Interactions")
@show ops
# --8<-- [end:read_opsum]



# --8<-- [start:read_permutation_group]
file = "triangular.9.hop.sublattices.tsl.toml"
fl = FileToml(file)
group = read_permutation_group(fl, "Symmetries")
@show group
# --8<-- [end:read_permutation_group]

# --8<-- [start:read_representation]
file = "irreps.toml"
fl = FileToml(file)

k_0 = read_representation(fl, "k_0")
@show k_0
@show isreal(k_0)

k_pi2 = read_representation(fl, "k_pi2")
@show k_pi2
@show isreal(k_pi2)

k_pi = read_representation(fl, "k_pi")
@show k_pi
@show isreal(k_pi)

k_pi2_half = read_representation(fl, "k_pi2_half")
@show k_pi2_half
@show isreal(k_pi2_half)
# --8<-- [end:read_representation]

