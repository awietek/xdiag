// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply_permutation.hpp"

#include <cstdint>

namespace xdiag::symmetries {

template <typename bit_t, typename int_t>
bit_t apply_permutation(bit_t state, int_t nsites, const int_t *permutation) {
  bit_t tstate = 0;
  for (int_t site = 0; site < nsites; ++site) {
    tstate |= ((state >> site) & 1) << permutation[site];
  }
  return tstate;
}

template uint16_t apply_permutation<uint16_t>(uint16_t, int32_t,
                                              const int64_t *);
template uint32_t apply_permutation<uint32_t>(uint32_t, int32_t,
                                              const int64_t *);
template uint64_t apply_permutation<uint64_t>(uint64_t, int32_t,
                                              const int64_t *);

template uint16_t apply_permutation<uint16_t>(uint16_t, int64_t,
                                              const int64_t *);
template uint32_t apply_permutation<uint32_t>(uint32_t, int64_t,
                                              const int64_t *);
template uint64_t apply_permutation<uint64_t>(uint64_t, int64_t,
                                              const int64_t *);
} // namespace xdiag::symmetries
