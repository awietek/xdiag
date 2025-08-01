// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "time_evolve_expokit.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/sparse/apply.hpp>
#include <xdiag/algebra/sparse/logic.hpp>
#include <xdiag/algorithms/time_evolution/zahexpv.hpp>
#include <xdiag/operators/logic/hc.hpp>
#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag {

template <typename op_t>
TimeEvolveExpokitResult
time_evolve_expokit(op_t const &ops, State state, double time, double precision,
                    int64_t m, double anorm, int64_t nnorm) try {
  auto res =
      time_evolve_expokit_inplace(ops, state, time, precision, m, anorm, nnorm);
  return {res.error, res.hump, state};
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

TimeEvolveExpokitResult time_evolve_expokit(OpSum const &ops, State state,
                                            double time, double precision,
                                            int64_t m, double anorm,
                                            int64_t nnorm) try {
  return time_evolve_expokit<OpSum>(ops, state, time, precision, m, anorm,
                                    nnorm);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename idx_t, typename coeff_t>
TimeEvolveExpokitResult
time_evolve_expokit(CSRMatrix<idx_t, coeff_t> const &ops, State state,
                    double time, double precision, int64_t m, double anorm,
                    int64_t nnorm) try {
  return time_evolve_expokit<CSRMatrix<idx_t, coeff_t>>(
      ops, state, time, precision, m, anorm, nnorm);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template TimeEvolveExpokitResult
time_evolve_expokit(CSRMatrix<int32_t, double> const &, State, double, double,
                    int64_t, double, int64_t);
template TimeEvolveExpokitResult
time_evolve_expokit(CSRMatrix<int32_t, complex> const &, State, double, double,
                    int64_t, double, int64_t);
template TimeEvolveExpokitResult
time_evolve_expokit(CSRMatrix<int64_t, double> const &, State, double, double,
                    int64_t, double, int64_t);
template TimeEvolveExpokitResult
time_evolve_expokit(CSRMatrix<int64_t, complex> const &, State, double, double,
                    int64_t, double, int64_t);

template <typename op_t>
TimeEvolveExpokitInplaceResult
time_evolve_expokit_inplace(op_t const &ops, State &state, double time,
                            double precision, int64_t m, double anorm,
                            int64_t nnorm) try {
  if (!ishermitian(ops)) {
    XDIAG_THROW("Input operator is not hermitian. Evolution using the expokit "
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
    return dot(block, v, w);
  };

  auto v0 = state.vectorC(0, false);
  auto t0 = rightnow();

  auto [err, hump] =
      zahexpv(time, apply_A, dot_f, v0, anorm, precision / time, m);

  timing(t0, rightnow(), "Time evolve expokit time", 1);
  return {err, hump};
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

TimeEvolveExpokitInplaceResult
time_evolve_expokit_inplace(OpSum const &ops, State &state, double time,
                            double precision, int64_t m, double anorm,
                            int64_t nnorm) try {
  return time_evolve_expokit_inplace<OpSum>(ops, state, time, precision, m,
                                            anorm, nnorm);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename idx_t, typename coeff_t>
TimeEvolveExpokitInplaceResult
time_evolve_expokit_inplace(CSRMatrix<idx_t, coeff_t> const &ops, State &state,
                            double time, double precision, int64_t m,
                            double anorm, int64_t nnorm) try {
  return time_evolve_expokit_inplace<CSRMatrix<idx_t, coeff_t>>(
      ops, state, time, precision, m, anorm, nnorm);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template TimeEvolveExpokitInplaceResult
time_evolve_expokit_inplace(CSRMatrix<int32_t, double> const &, State &, double,
                            double, int64_t, double, int64_t);
template TimeEvolveExpokitInplaceResult
time_evolve_expokit_inplace(CSRMatrix<int32_t, complex> const &, State &,
                            double, double, int64_t, double, int64_t);
template TimeEvolveExpokitInplaceResult
time_evolve_expokit_inplace(CSRMatrix<int64_t, double> const &, State &, double,
                            double, int64_t, double, int64_t);
template TimeEvolveExpokitInplaceResult
time_evolve_expokit_inplace(CSRMatrix<int64_t, complex> const &, State &,
                            double, double, int64_t, double, int64_t);

} // namespace xdiag
