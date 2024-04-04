#include "time_evolution.h"

#include <xdiag/algebra/algebra.h>
#include <xdiag/algebra/apply.h>

#include <xdiag/algorithms/norm_estimate.h>
#include <xdiag/algorithms/time_evolution/exp_sym_v.h>
#include <xdiag/algorithms/time_evolution/zahexpv.h>

#include <xdiag/utils/timing.h>

namespace xdiag {

std::tuple<double, double> time_evolve_inplace(BondList const &bonds,
                                               State &state, double time,
                                               double precision, int64_t m,
                                               double anorm,
                                               int64_t nnorm) try {
  if (state.isreal()) {
    state.make_complex();
  }
  auto const &block = state.block();

  if (anorm == 0.) { // if anorm is default value 0., compute an estimate
    for (int64_t j = 0; j < nnorm; ++j) {
      double anormj = norm_estimate(bonds, block);
      if (anormj > anorm) {
        anorm = anormj;
      }
    }
    Log(1, "norm estimate: {}", anorm);
  }

  int64_t iter = 1;
  auto apply_A = [&iter, &bonds, &block](arma::cx_vec const &v) {
    auto ta = rightnow();
    auto w = arma::cx_vec(v.n_rows, arma::fill::zeros);
    apply(bonds, block, v, block, w);
    w *= -1.0i;
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

  timing(t0, rightnow(), "Zahexpv time", 1);
  return {err, hump};
} catch (...) {
  XDiagRethrow("Unable to perform inplace time evolution");
  return {0., 0.};
}

State time_evolve(BondList const &bonds, State state, double time,
                  double precision, int64_t m, double anorm,
                  int64_t nnorm) try {
  auto [err, hump] =
      time_evolve_inplace(bonds, state, time, precision, m, anorm, nnorm);
  Log(2, "error (estimated): {}, hump: {}", err, hump);
  return state;
} catch (...) {
  XDiagRethrow("Unable to perform real time evolution");
  return State();
}

State imag_time_evolve(BondList const &bonds, State const &state, double time,
                       double precision = 1e-12, int64_t max_iterations = 1000,
                       double deflation_tol = 1e-7) try {
  return exp_sym_v(bonds, state, -time, false, true, precision, max_iterations,
                   deflation_tol);
} catch (...) {
  XDiagRethrow("Unable to perform imaginary time evolution");
  return State();
}

} // namespace xdiag
