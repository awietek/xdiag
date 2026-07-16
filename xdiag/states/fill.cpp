// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fill.hpp"

#include <variant>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/random/hash.hpp>
#include <xdiag/random/hash_functions.hpp>
#include <xdiag/random/random_utils.hpp>
#include <xdiag/states/norm.hpp>
#include <xdiag/utils/error.hpp>

#ifdef XDIAG_DISTRIBUTED
#include <mpi.h>
#endif

namespace xdiag {

template <typename coeff_t, class block_t, class coeff_f>
void fill(block_t const &block, arma::Col<coeff_t> &vec, coeff_f coeff) try {
  int64_t idx = 0;
  for (auto const &pstate : block) {
    vec(idx++) = coeff(pstate);
  }
}
XDIAG_CATCH

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
}
XDIAG_CATCH

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
}
XDIAG_CATCH

void fill(State &state, RandomState const &rstate, int64_t col) try {
  int64_t seed = rstate.seed();
  int64_t seed_modified =
      random::hash_combine(seed, random::hash(state.block()));

#ifdef XDIAG_DISTRIBUTED
  // For distributed blocks every rank holds a disjoint slice of the global
  // vector, but fill_random_normal_vector always starts its (seeded) random
  // sequence at local index 0. Seeding all ranks identically would therefore
  // fill those slices with the *same* local sequence, producing a global start
  // vector that is correlated across ranks -- and for tiny blocks (e.g. one
  // amplitude per rank) an exact symmetry eigenvector, onto which Lanczos then
  // locks (returning the wrong eigenvalue). Offsetting the seed by the MPI rank
  // makes each rank's slice an independent draw, so the global vector is a
  // genuine random vector.
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if ((mpi_size > 1) && isdistributed(state.block())) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    seed_modified = random::hash_combine((uint64_t)seed_modified,
                                         random::hash_fnv1((uint64_t)rank));
  }
#endif

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
}
XDIAG_CATCH

// Set the given column to the basis vector of a single product state: 1 on the
// matching basis element, 0 elsewhere. Implemented by reusing the coefficient
// fill, which iterates the block's product states, so it is correct for every
// basis type (including symmetric bases, where the match is against the
// canonical representative product state).
void fill(State &state, ProductState const &pstate, int64_t col) try {
  if (state.nsites() != pstate.size()) {
    XDIAG_THROW("State and ProductState do not have the same number of sites");
  } else if (col >= state.ncols()) {
    XDIAG_THROW("Column index larger than number of columns in State");
  }
  std::function<double(ProductState const &)> coeff =
      [&pstate](ProductState const &p) -> double {
    return (p == pstate) ? 1.0 : 0.0;
  };
  fill(state, coeff, col);
}
XDIAG_CATCH

// void fill(State &state, GPWF const &gpwf, int64_t col) try {
//   if (state.nsites() != gpwf.nsites()) {
//     XDIAG_THROW("State and GPWF do not have the same number of sites");
//   }
//   auto const &block = state.block();

// #ifdef XDIAG_DISTRIBUTED
//   if (!std::holds_alternative<Spinhalf>(block) &&
//       !std::holds_alternative<SpinhalfDistributed>(block)) {
//     XDIAG_THROW("GPWF is currently only defined for \"Spinhalf\" and "
//                 "\"SpinhalfDistributed\" type blocks");
//   }
// #else
//   if (!std::holds_alternative<Spinhalf>(block)) {
//     XDIAG_THROW("GPWF is currently only defined for \"Spinhalf\" type blocks");
//   }
// #endif

//   if (gpwf.isreal()) {
//     std::function<double(ProductState const &)> f =
//         [&](ProductState const &pstate) { return gpwf.coefficient(pstate); };
//     fill(state, f, col);
//   } else {
//     std::function<complex(ProductState const &)> f =
//         [&](ProductState const &pstate) { return gpwf.coefficientC(pstate); };
//     state.make_complex();
//     fill(state, f, col);
//   }
// }
// XDIAG_CATCH

} // namespace xdiag
