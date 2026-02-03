// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "is_valid_permutation.hpp"

#include <algorithm>
#include <cstdint>

namespace xdiag::symmetries {

template <typename int_t>
bool is_valid_permutation(int_t nsites, const int_t *permutation) {
  for (int_t i = 0; i < nsites; ++i) {
    if (std::find(permutation, permutation + nsites, i) == permutation + nsites)
      return false;
  }
  return true;
}

template bool is_valid_permutation(int32_t nsites, const int32_t *permutation);
template bool is_valid_permutation(int64_t nsites, const int64_t *permutation);

} // namespace xdiag::symmetries
