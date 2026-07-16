// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "time_evolve_expokit.hpp"

#include <xdiag/algebra/ishermitian.hpp>
#include <xdiag/kernels/apply.hpp>
#include <xdiag/kernels/sparse/apply.hpp>
#include <xdiag/linalg/time_evolution/zahexpv.hpp>
#include <xdiag/math/dot.hpp>
#include <xdiag/math/norm.hpp>
#include <xdiag/states/norm.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag {

template <typename op_t>
TimeEvolveExpokitResult
time_evolve_expokit(op_t const &ops, State state, double time, double precision,
                    int64_t m, double anorm, int64_t nnorm) try {
  auto res =
      time_evolve_expokit_inplace(ops, state, time, precision, m, anorm, nnorm);
  return {res.error, res.hump, state};
}
XDIAG_CATCH

TimeEvolveExpokitResult time_evolve_expokit(OpSum const &ops, State state,
                                            double time, double precision,
                                            int64_t m, double anorm,
                                            int64_t nnorm) try {
  return time_evolve_expokit<OpSum>(ops, state, time, precision, m, anorm,
                                    nnorm);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
TimeEvolveExpokitResult
time_evolve_expokit(CSRMatrix<idx_t, coeff_t> const &ops, State state,
                    double time, double precision, int64_t m, double anorm,
                    int64_t nnorm) try {
  return time_evolve_expokit<CSRMatrix<idx_t, coeff_t>>(
      ops, state, time, precision, m, anorm, nnorm);
}
XDIAG_CATCH

#define XDIAG_INST(IDX, COEFF)                                                 \
  template TimeEvolveExpokitResult time_evolve_expokit(                        \
      CSRMatrix<IDX, COEFF> const &, State, double, double, int64_t, double,   \
      int64_t);
XDIAG_INST(int32_t, double)
XDIAG_INST(int32_t, complex)
XDIAG_INST(int64_t, double)
XDIAG_INST(int64_t, complex)
#undef XDIAG_INST

template <typename op_t>
TimeEvolveExpokitInplaceResult
time_evolve_expokit_inplace(op_t const &ops, State &state, double time,
                            double precision, int64_t m, double anorm,
                            int64_t nnorm) try {

  if (dim(state) == 0) {
    Log.warn("Warning: initial state zero dimensional in "
             "time_evolve_expokit_inplace");
    return TimeEvolveExpokitInplaceResult();
  }

  if (!ishermitian(ops, state.block())) {
    XDIAG_THROW("Input OpSum is not hermitian. Evolution using the expokit "
                "algorithm requires the operator to be hermitian.");
  }
  if (!isvalid(state)) {
    XDIAG_THROW("Initial state must be a valid state (i.e. not default "
                "constructed by e.g. an annihilation operator)");
  }

  if (norm(state) == 0.) {
    XDIAG_THROW("Initial state has zero norm");
  }

  if (state.isreal()) {
    state.make_complex();
  }
  auto const &block = state.block();

  if (anorm == 0.) { // if anorm is default value 0., compute an estimate
    for (int64_t j = 0; j < nnorm; ++j) {
      double anormj = norm_estimate(ops, block);
      if (anormj > anorm) {
        anorm = anormj;
      }
    }
    Log(1, "norm estimate: {}", anorm);
  }

  int64_t iter = 1;
  auto apply_A = [&iter, &ops, &block](arma::cx_vec const &v) {
    auto ta = rightnow();
    auto w = arma::cx_vec(v.n_rows, arma::fill::zeros);
    apply(ops, block, v, block, w);
    w *= complex(0.0, -1.0);
    Log(2, "Lanczos iteration {}", iter);
    timing(ta, rightnow(), "MVM", 2);
    ++iter;
    return w;
  };
  auto dot_f = [&block](arma::cx_vec const &v, arma::cx_vec const &w) {
    return math::dot(block, v, w);
  };

  auto v0 = state.vectorC(0, false);
  auto t0 = rightnow();

  auto [err, hump] =
      zahexpv(time, apply_A, dot_f, v0, anorm, precision / time, m);

  timing(t0, rightnow(), "Time evolve expokit time", 1);
  return {err, hump};
}
XDIAG_CATCH

TimeEvolveExpokitInplaceResult
time_evolve_expokit_inplace(OpSum const &ops, State &state, double time,
                            double precision, int64_t m, double anorm,
                            int64_t nnorm) try {
  return time_evolve_expokit_inplace<OpSum>(ops, state, time, precision, m,
                                            anorm, nnorm);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
TimeEvolveExpokitInplaceResult
time_evolve_expokit_inplace(CSRMatrix<idx_t, coeff_t> const &ops, State &state,
                            double time, double precision, int64_t m,
                            double anorm, int64_t nnorm) try {
  return time_evolve_expokit_inplace<CSRMatrix<idx_t, coeff_t>>(
      ops, state, time, precision, m, anorm, nnorm);
}
XDIAG_CATCH

#define XDIAG_INST(IDX, COEFF)                                                 \
  template TimeEvolveExpokitInplaceResult time_evolve_expokit_inplace(         \
      CSRMatrix<IDX, COEFF> const &, State &, double, double, int64_t, double, \
      int64_t);
XDIAG_INST(int32_t, double)
XDIAG_INST(int32_t, complex)
XDIAG_INST(int64_t, double)
XDIAG_INST(int64_t, complex)
#undef XDIAG_INST

} // namespace xdiag
