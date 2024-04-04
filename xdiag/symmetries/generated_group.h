#pragma once

#include <utility>

#include <xdiag/symmetries/permutation.h>
#include <xdiag/symmetries/permutation_group.h>
#include <xdiag/symmetries/representation.h>

namespace xdiag {

PermutationGroup generated_group(Permutation const &generator);
PermutationGroup generated_group(std::vector<Permutation> const &generators);

Representation generated_irrep(Permutation const &generator, complex phase);
Representation generated_irrep(std::vector<Permutation> const &generators,
                               std::vector<complex> const &phases);

std::pair<PermutationGroup, Representation>
generated_group_irrep(Permutation const &generator, complex phase);
std::pair<PermutationGroup, Representation>
generated_group_irrep(std::vector<Permutation> const &generators,
                      std::vector<complex> const &phases);

} // namespace xdiag
