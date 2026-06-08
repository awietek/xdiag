// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "valid.hpp"

#include <cstdint>
#include <limits>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag {

template <typename idx_t, typename coeff_t, typename block_t>
void check_valid_sparse_matrix(OpSum const &ops, block_t const &block_in,
                               block_t block_out, idx_t i0) try {

  if (!((i0 == 0) || (i0 == 1))) {
    XDIAG_THROW(fmt::format(
        "Invalid zero index i0. Must be either 0 or 1, but got i0={}", i0));
  }

  // Check if ops and blocks are compatible
  if (!blocks_match(ops, block_in, block_out)) {
    XDIAG_THROW("Cannot create a sparse matrix on blocks. The resulting block "
                "is not in "
                "the correct symmetry sector. Please check the quantum numbers "
                "of the output block.");
  }

  // Check if real matrix can be created
  if constexpr (isreal<coeff_t>()) {
    if (!isreal(ops)) {
      XDIAG_THROW(
          "Cannot create a real sparse matrix from an Op or OpSum which "
          "is complex.");
    }
    if (!isreal(block_in) || !isreal(block_out)) {
      XDIAG_THROW("Cannot create a real sparse matrix when a block is complex.")
    }
  }

  int64_t m = block_out.size();
  int64_t n = block_in.size();
  if ((m > std::numeric_limits<idx_t>::max()) ||
      (n > std::numeric_limits<idx_t>::max())) {
    XDIAG_THROW(
        "Block is too large for index type which attempts to hold indices. "
        "Consider using a larger index type, e.g. 64 bit integers");
  }
}
XDIAG_CATCH

// Explicit instantiations for each (idx_t, coeff_t) over the supported blocks.
#define XDIAG_INSTANTIATE_CHECK_VALID(idx_t, coeff_t, block_t)                 \
  template void check_valid_sparse_matrix<idx_t, coeff_t, block_t>(            \
      OpSum const &, block_t const &, block_t, idx_t);

#define XDIAG_INSTANTIATE_CHECK_VALID_BLOCK(block_t)                           \
  XDIAG_INSTANTIATE_CHECK_VALID(int32_t, double, block_t)                      \
  XDIAG_INSTANTIATE_CHECK_VALID(int32_t, complex, block_t)                     \
  XDIAG_INSTANTIATE_CHECK_VALID(int64_t, double, block_t)                      \
  XDIAG_INSTANTIATE_CHECK_VALID(int64_t, complex, block_t)

XDIAG_INSTANTIATE_CHECK_VALID_BLOCK(Spinhalf)
XDIAG_INSTANTIATE_CHECK_VALID_BLOCK(Boson)

#undef XDIAG_INSTANTIATE_CHECK_VALID_BLOCK
#undef XDIAG_INSTANTIATE_CHECK_VALID

} // namespace xdiag
