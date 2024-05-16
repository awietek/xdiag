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
