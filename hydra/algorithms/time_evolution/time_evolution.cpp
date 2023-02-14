#include "time_evolution.h"

#include <hydra/algebra/algebra.h>
#include <hydra/algorithms/time_evolution/exp_sym_v.h>
#include <hydra/algorithms/time_evolution/zahexpv.h>
#include <hydra/utils/timing.h>

namespace hydra {

std::tuple<double, double> time_evolve_inplace(BondList const &bonds,
                                               StateCplx &state, double time,
                                               double precision, int m,
                                               double anorm, int nnorm) {
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

  auto &v0 = state.vector();
  auto t0 = rightnow();
  auto [err, hump] =
      zahexpv(time, apply_A, v0, precision / time, m, anorm, nnorm);
  timing(t0, rightnow(), "Zahexpv time", 1);
  return {err, hump};
}

template <>
StateCplx time_evolve(BondList const &bonds, StateReal state, double time,
                      double precision, int m, double anorm, int nnorm) {
  auto state_cplx = to_cplx(state);
  auto [err, hump] = time_evolve_inplace(bonds, state_cplx, time, precision);
  Log(2, "error (estimated): {}, hump: {}", err, hump);
  return state_cplx;
}

template <>
StateCplx time_evolve(BondList const &bonds, StateCplx state, double time,
                      double precision, int m, double anorm, int nnorm) {
  auto [err, hump] = time_evolve_inplace(bonds, state, time, precision);
  Log(2, "error (estimated): {}, hump: {}", err, hump);
  return state;
}

template <typename coeff_t>
double imag_time_evolve_inplace(BondList const &bonds, State<coeff_t> &state,
                                double time, double precision) {
  return exp_sym_v_inplace(bonds, state, -time, precision);
}

template double imag_time_evolve_inplace(BondList const &, State<double> &,
                                         double, double);
template double imag_time_evolve_inplace(BondList const &, State<complex> &,
                                         double, double);

template <typename coeff_t>
State<coeff_t> imag_time_evolve(BondList const &bonds, State<coeff_t> state,
                                double time, double precision) {
  imag_time_evolve_inplace(bonds, state, time, precision);
  return state;
}
template State<double> imag_time_evolve(BondList const &, State<double>, double,
                                        double);
template State<complex> imag_time_evolve(BondList const &, State<complex>,
                                         double, double);

} // namespace hydra
