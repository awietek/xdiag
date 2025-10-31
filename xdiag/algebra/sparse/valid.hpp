// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>
#include <xdiag/operators/logic/block.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

template <typename idx_t, typename coeff_t, typename block_t>
inline void check_valid_sparse_matrix(OpSum const &ops, block_t const &block_in,
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

  int64_t nsites = block_in.nsites();
  check_valid(ops, nsites);

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
} // namespace xdiag
