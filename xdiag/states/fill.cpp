// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fill.hpp"

#include <variant>

#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/random/hash.hpp>
#include <xdiag/random/hash_functions.hpp>
#include <xdiag/random/random_utils.hpp>

namespace xdiag {

template <typename coeff_t, class block_t, class coeff_f>
void fill(block_t const &block, arma::Col<coeff_t> &vec, coeff_f coeff) try {
  int64_t idx = 0;
  for (auto const &pstate : block) {
    vec(idx++) = coeff(pstate);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void fill(State &state, std::function<double(ProductState const &)> coeff_f,
          int64_t col) try {
  auto const &block = state.block();
  if (state.isreal()) {
    arma::vec v = state.vector(col, false);
    std::visit([&](auto &&block) { fill(block, v, coeff_f); }, block);
  } else {
    arma::cx_vec v = state.vectorC(col, false);
    auto coeff_f_c = [&](ProductState const &pstate) {
      return (complex)coeff_f(pstate);
    };
    std::visit([&](auto &&block) { fill(block, v, coeff_f_c); }, block);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void fill(State &state, std::function<complex(ProductState const &)> coeff_f,
          int64_t col) try {
  auto const &block = state.block();
  if (state.isreal()) {
    XDIAG_THROW("Cannot fill real state with complex coefficients. Maybe use "
                "\"make_complex\" first?");
  } else {
    arma::cx_vec v = state.vectorC(col, false);
    std::visit([&](auto &&block) { fill(block, v, coeff_f); }, block);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void fill(State &state, RandomState const &rstate, int64_t col) try {
  int64_t seed = rstate.seed();
  int64_t seed_modified =
      random::hash_combine(seed, random::hash(state.block()));

  if (state.isreal()) {
    auto v = state.vector(col, false);
    random::fill_random_normal_vector(v, seed_modified);
  } else {
    auto v = state.vectorC(col, false);
    random::fill_random_normal_vector(v, seed_modified);
  }
  if (rstate.normalized()) {
    double nrm = norm(state);
    state /= nrm;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <class block_t, typename coeff_t>
void fill(block_t const &block, arma::Col<coeff_t> &vec,
          ProductState const &p) try {
  int64_t idx = block.index(p);
  if (isdistributed(block)) {
    vec.zeros();

#ifdef XDIAG_USE_MPI
    // all processes except one should hold invalid index
    int64_t max_idx;
    MPI_Allreduce(&idx, &max_idx, 1, MPI_LONG_LONG_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    if (max_idx == invalid_index) {
      XDIAG_THROW("Index of product state cannot be determined");
    } else if (idx != invalid_index) {
      vec(idx) = 1.0;
    }
#endif
  } else {
    if (idx == invalid_index) {
      XDIAG_THROW("Index of product state cannot be determined");
    }
    vec.zeros();
    vec(idx) = 1.0;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void fill(State &state, ProductState const &pstate, int64_t col) try {
  if (state.nsites() != pstate.size()) {
    XDIAG_THROW("State and ProductState do not have the same number of sites");
  } else if (col >= state.ncols()) {
    XDIAG_THROW("Column index larger than number of columns in State");
  }
  auto const &block = state.block();
  if (state.isreal()) {
    arma::vec v = state.vector(col, false);
    std::visit([&](auto &&block) { fill(block, v, pstate); }, block);
  } else {
    arma::cx_vec v = state.vectorC(col, false);
    std::visit([&](auto &&block) { fill(block, v, pstate); }, block);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void fill(State &state, GPWF const &gpwf, int64_t col) try {
  if (state.nsites() != gpwf.nsites()) {
    XDIAG_THROW("State and GPWF do not have the same number of sites");
  }
  auto const &block = state.block();

#ifdef XDIAG_USE_MPI
  if (!std::holds_alternative<Spinhalf>(block) &&
      !std::holds_alternative<SpinhalfDistributed>(block)) {
    XDIAG_THROW("GPWF is currently only defined for \"Spinhalf\" and "
                "\"SpinhalfDistributed\" type blocks");
  }
#else
  if (!std::holds_alternative<Spinhalf>(block)) {
    XDIAG_THROW("GPWF is currently only defined for \"Spinhalf\" type blocks");
  }
#endif

  if (gpwf.isreal()) {
    std::function<double(ProductState const &)> f =
        [&](ProductState const &pstate) { return gpwf.coefficient(pstate); };
    fill(state, f, col);
  } else {
    std::function<complex(ProductState const &)> f =
        [&](ProductState const &pstate) { return gpwf.coefficientC(pstate); };
    state.make_complex();
    fill(state, f, col);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
