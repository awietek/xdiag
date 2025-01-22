#pragma once
#include <string>
#include <utility>
#include <vector>

#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/extern/gsl/span>
#include <xdiag/symmetries/operations/fermi_sign.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::symmetries {

using span_size_t = gsl::span<int const>::size_type;

bool is_valid_permutation(int nsites, const int *permutation);

template <typename bit_t>
bit_t apply_permutation(bit_t state, int nsites, const int *permutation);

template <typename bit_t>
bit_t apply_permutation(bit_t state, gsl::span<int const> permutation);

} // namespace xdiag::symmetries
