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

void fill(State &state, ProductState const &pstate, int64_t col) try {
  if (state.n_sites() != pstate.n_sites()) {
    XDIAG_THROW("State and ProductState do not have the same number of sites");
  } else if (col >= state.n_cols()) {
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

template void fill(Spinhalf const &, arma::vec &, ProductState const &);
template void fill(Spinhalf const &, arma::cx_vec &, ProductState const &);
template void fill(tJ const &, arma::vec &, ProductState const &);
template void fill(tJ const &, arma::cx_vec &, ProductState const &);
template void fill(Electron const &, arma::vec &, ProductState const &);
template void fill(Electron const &, arma::cx_vec &, ProductState const &);
#ifdef XDIAG_USE_MPI
template void fill(SpinhalfDistributed const &, arma::vec &,
                   ProductState const &);
template void fill(SpinhalfDistributed const &, arma::cx_vec &,
                   ProductState const &);
template void fill(tJDistributed const &, arma::vec &, ProductState const &);
template void fill(tJDistributed const &, arma::cx_vec &, ProductState const &);
#endif


} // namespace xdiag
