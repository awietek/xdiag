// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "create_state.hpp"

#include <xdiag/math/dot.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {
template <typename block_t>
static State product_state(block_t const &block,
                           std::vector<int64_t> const &local_state,
                           bool real) try {
  State state(block, real);
  ProductState pstate(local_state);
  fill(state, pstate);
  return state;
}
XDIAG_CATCH

State product_state(Block const &block, std::vector<int64_t> const &local_state,
                    bool real) try {
  return std::visit(
      [&](auto &&block) { return product_state(block, local_state, real); },
      block);
}
XDIAG_CATCH

template <typename block_t>
static State random_state(block_t const &block, bool real, int64_t ncols,
                          int64_t seed, bool normalized) try {
  auto state = State(block, real, ncols);
  if (ncols == 1) {
    auto rstate = RandomState(seed, normalized);
    fill(state, rstate);
  } else {

    std::mt19937 rng(seed);
    std::uniform_int_distribution<std::mt19937::result_type> dist(0);

    // Fill all random columns
    for (int64_t col = 0; col < ncols; ++col) {
      auto rstate = RandomState(seed, false);
      fill(state, rstate, col);
      seed += dist(rng);
    }

    // orthonormalize
    if (normalized) {
      if (ncols > dim(block)) {
        XDIAG_THROW(
            "Cannot orthonormalize random state, since the block dimension is "
            "smaller than the requesed number of columns.  Either don't try to "
            "normalize the state or choose a smaller number of colums");
      }
      if (real) {
        auto V = state.matrix(false);
        auto VdagV = math::matrix_dot(block, V, V);
        if (VdagV.is_sympd()) {
          arma::mat L = chol(VdagV, "lower");
          V = V * inv(trimatl(L)).t();
        } else {
          XDIAG_THROW(
              "Random vectors cannot be orthonormalized. The appear to be "
              "linearly dependent (please report, this should not happen)");
        }
      } else { // complex
        auto V = state.matrixC(false);
        auto VdagV = math::matrix_dot(block, V, V);
        if (VdagV.is_sympd()) {
          arma::cx_mat L = chol(VdagV, "lower");
          V = V * inv(trimatl(L)).t();
        } else {
          XDIAG_THROW(
              "Random vectors cannot be orthonormalized. The appear to be "
              "linearly dependent (please report, this should not happen)");
        }
      }
    }
  }
  return state;
}
XDIAG_CATCH

State random_state(Block const &block, bool real, int64_t ncols, int64_t seed,
                   bool normalized) try {
  return std::visit(
      [&](auto &&block) {
        return random_state(block, real, ncols, seed, normalized);
      },
      block);
}
XDIAG_CATCH

template <typename block_t>
static State zero_state(block_t const &block, bool real, int64_t ncols) {
  return State(block, real, ncols);
}

State zero_state(Block const &block, bool real, int64_t ncols) {
  return std::visit([&](auto &&blk) { return zero_state(blk, real, ncols); },
                    block);
}

void zero(State &state) {
  if (state.isreal()) {
    state.matrix(false).zeros();
  } else {
    state.matrixC(false).zeros();
  }
}

} // namespace xdiag
