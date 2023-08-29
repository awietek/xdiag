#include "time_evolution.h"

#include <hydra/algebra/algebra.h>
#include <hydra/algebra/apply.h>

#include <hydra/algorithms/time_evolution/exp_sym_v.h>
#include <hydra/algorithms/time_evolution/zahexpv.h>
#include <hydra/utils/timing.h>

namespace hydra {

std::tuple<double, double> time_evolve_inplace(BondList const &bonds,
                                               State &state, double time,
                                               double precision, int m,
                                               double anorm, int nnorm) {
  if (state.isreal()) {
    state.make_complex();
  }
  auto const &block = state.block();

  int iter = 1;

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

  auto v0 = state.vectorC(0, false);
  auto t0 = rightnow();
  auto [err, hump] =
      zahexpv(time, apply_A, v0, precision / time, m, anorm, nnorm);
  timing(t0, rightnow(), "Zahexpv time", 1);
  return {err, hump};
}

State time_evolve(BondList const &bonds, State state, double time,
                  double precision, int m, double anorm, int nnorm) {
  auto [err, hump] =
      time_evolve_inplace(bonds, state, time, precision, m, anorm, nnorm);
  Log(2, "error (estimated): {}, hump: {}", err, hump);
  return state;
}

State imag_time_evolve(BondList const &bonds, State const &state, double time,
                       double precision = 1e-12, int64_t max_iterations = 1000,
                       double deflation_tol = 1e-7) {
  return exp_sym_v(bonds, state, -time, false, true, precision, max_iterations,
                   deflation_tol);
}

} // namespace hydra
