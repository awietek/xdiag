// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "create_state.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/utils/xdiag_show.hpp>

namespace xdiag {
State product_state(Block const &block,
                    std::vector<std::string> const &local_state,
                    bool real) try {
  return std::visit(
      [&](auto &&block) { return product_state(block, local_state, real); },
      block);
}
XDIAG_CATCH

template <typename block_t>
State product_state(block_t const &block,
                    std::vector<std::string> const &local_state,
                    bool real) try {
  auto state = State(block, real);
  auto pstate = ProductState(local_state);
  fill(state, pstate);
  return state;
}
XDIAG_CATCH

template State product_state(Spinhalf const &, std::vector<std::string> const &,
                             bool);
template State product_state(tJ const &, std::vector<std::string> const &,
                             bool);
template State product_state(Electron const &, std::vector<std::string> const &,
                             bool);

#ifdef XDIAG_USE_MPI
template State product_state(SpinhalfDistributed const &,
                             std::vector<std::string> const &, bool);
template State product_state(tJDistributed const &,
                             std::vector<std::string> const &, bool);
template State product_state(ElectronDistributed const &,
                             std::vector<std::string> const &, bool);
#endif

State random_state(Block const &block, bool real, int64_t ncols, int64_t seed,
                   bool normalized) try {
  Log("asdf");
  return std::visit(
      [&](auto &&block) {
        return random_state(block, real, ncols, seed, normalized);
      },
      block);
}
XDIAG_CATCH

template <typename block_t>
State random_state(block_t const &block, bool real, int64_t ncols, int64_t seed,
                   bool normalized) try {
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
      if (ncols > size(block)) {
        XDIAG_THROW(
            "Cannot orthonormalize random state, since the block dimension is "
            "smaller than the requesed number of columns.  Either don't try to "
            "normalize the state or choose a smaller number of colums");
      }
      if (real) {
        auto V = state.matrix(false);
        auto VdagV = matrix_dot(block, V, V);
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
        auto VdagV = matrix_dot(block, V, V);
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

template State random_state(Spinhalf const &, bool, int64_t, int64_t, bool);
template State random_state(tJ const &, bool, int64_t, int64_t, bool);
template State random_state(Electron const &, bool, int64_t, int64_t, bool);
#ifdef XDIAG_USE_MPI
template State random_state(SpinhalfDistributed const &, bool, int64_t, int64_t,
                            bool);
template State random_state(tJDistributed const &, bool, int64_t, int64_t,
                            bool);
template State random_state(ElectronDistributed const &, bool, int64_t, int64_t,
                            bool);
#endif

State zero_state(Block const &block, bool real, int64_t ncols) {
  return std::visit([&](auto &&blk) { return zero_state(blk, real, ncols); },
                    block);
}

template <typename block_t>
State zero_state(block_t const &block, bool real, int64_t ncols) {
  return State(block, real, ncols);
}
template State zero_state(Spinhalf const &, bool, int64_t);
template State zero_state(tJ const &, bool, int64_t);
template State zero_state(Electron const &, bool, int64_t);
#ifdef XDIAG_USE_MPI
template State zero_state(ElectronDistributed const &, bool, int64_t);
template State zero_state(tJDistributed const &, bool, int64_t);
template State zero_state(SpinhalfDistributed const &, bool, int64_t);
#endif

void zero(State &state) {
  if (state.isreal()) {
    state.matrix(false).zeros();
  } else {
    state.matrixC(false).zeros();
  }
}

} // namespace xdiag
