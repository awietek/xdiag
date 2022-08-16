#pragma once

#include <lila/external/gsl/span>

#include <hydra/combinatorics/combinations.h>
#include <hydra/symmetries/operations/fermi_sign.h>
#include <hydra/symmetries/representation.h>
#include <hydra/symmetries/permutation.h>

#include <string>
#include <utility>
#include <vector>

namespace hydra {
std::vector<Permutation> read_permutations(std::string filename);
}

namespace hydra::symmetries {

using span_size_t = gsl::span<int const>::size_type;

bool is_valid_permutation(int n_sites, const int *permutation);

template <typename bit_t>
bit_t apply_permutation(bit_t state, int n_sites, const int *permutation);

template <typename bit_t>
bit_t apply_permutation(bit_t state, gsl::span<int const> permutation);

} // namespace hydra::symmetries
