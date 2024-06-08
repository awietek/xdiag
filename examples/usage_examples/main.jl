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
irrep = Representation([1, -1, 1, -1])
block_sym = Spinhalf(N, group, irrep)
@show block_sym

# with symmetries and Sz
block_sym_sz = Spinhalf(N, nup, group, irrep)
@show block_sym_sz

@show n_sites(block_sym_sz)
@show size(block_sym_sz)
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
irrep = Representation([1, -1, 1, -1])
block_sym = tJ(N, nup, ndn, group, irrep)
@show block_sym

@show n_sites(block_sym)
@show size(block_sym)
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
irrep = Representation([1, -1, 1, -1])
block_sym = Electron(N, group, irrep)
@show block_sym

# with symmetries and number conservation
block_sym_np = Electron(N, nup, ndn, group, irrep)
@show block_sym_np

@show n_sites(block_sym_np)
@show size(block_sym_np)
# --8<-- [end:Electron]


# --8<-- [start:matrix]
let
    # Creates matrix H_{k=2} in Eq (18.23) of https://link.springer.com/content/pdf/10.1007/978-3-540-74686-7_18.pdf
    N = 4
    nup = 3
    ndn = 2

    # Define a Hubbard chain model
    bonds = BondList()
    for i in 1:N
        bonds += Bond("HOP", "T", [i-1, i % N])
    end
    bonds["T"] = 1.0;
    bonds["U"] = 5.0;

    # Create the a permutation group
    p1 = Permutation([0, 1, 2, 3])
    p2 = Permutation([1, 2, 3, 0])
    p3 = Permutation([2, 3, 0, 1])
    p4 = Permutation([3, 0, 1, 2])
    group = PermutationGroup([p1, p2, p3, p4])
    irrep = Representation([1, -1, 1, -1])
    block = Electron(N, nup, ndn, group, irrep)

    H = matrix(bonds, block)
    display(H)
end
# --8<-- [end:matrix]

# --8<-- [start:eigval0]
let 
    N = 8;
    nup = N ÷ 2;
    block = Spinhalf(N, nup);
    
    # Define the nearest-neighbor Heisenberg model
    bonds = BondList()
    for i in 1:N
        bonds += Bond("HB", "J", [i-1, i % N])
    end
    bonds["J"] = 1.0;

    e0 = eigval0(bonds, block);
end
# --8<-- [end:eigval0]

# --8<-- [start:eig0]
let 
    N = 8;
    nup = N ÷ 2;
    block = Spinhalf(N, nup);
    
    # Define the nearest-neighbor Heisenberg model
    bonds = BondList()
    for i in 1:N
        bonds += Bond("HB", "J", [i-1, i % N])
    end
    bonds["J"] = 1.0;

    e0, gs = eig0(bonds, block);
end
# --8<-- [end:eig0]
