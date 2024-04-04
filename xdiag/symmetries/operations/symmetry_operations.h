#pragma once

#include <xdiag/extern/gsl/span>

#include <xdiag/combinatorics/combinations.h>
#include <xdiag/symmetries/operations/fermi_sign.h>
#include <xdiag/symmetries/representation.h>
#include <xdiag/symmetries/permutation.h>

#include <string>
#include <utility>
#include <vector>

namespace xdiag {
std::vector<Permutation> read_permutations(std::string filename);
}

namespace xdiag::symmetries {

using span_size_t = gsl::span<int const>::size_type;

bool is_valid_permutation(int n_sites, const int *permutation);

template <typename bit_t>
bit_t apply_permutation(bit_t state, int n_sites, const int *permutation);

template <typename bit_t>
bit_t apply_permutation(bit_t state, gsl::span<int const> permutation);

} // namespace xdiag::symmetries
