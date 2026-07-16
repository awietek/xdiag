# SPDX-License-Identifier: Apache-2.0
# STATIC hand-written special -- copied verbatim by generate.sh.
#
# PermutationGroup from a Vector{Permutation}: the C++ vector<Permutation>
# constructor cannot be marshalled by the generator, so we assemble the
# (nsites x n_permutations) matrix the generated Matrix{Int64} constructor
# expects (0-based, one permutation per column) and forward to it.

function PermutationGroup(perms::Vector{Permutation})
    isempty(perms) && return PermutationGroup()
    nsites = size(perms[1])
    mat = Matrix{Int64}(undef, nsites, length(perms))   # column j = permutation j
    for (j, perm) in enumerate(perms)
        for i in 1:nsites
            mat[i, j] = perm[i] - 1                      # perm[i] is 1-based
        end
    end
    return PermutationGroup(mat)
end


