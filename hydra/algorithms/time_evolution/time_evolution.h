#pragma once

#include <hydra/algorithms/time_evolution/exp_sym_v.h>
#include <hydra/algorithms/time_evolution/zahexpv.h>
#include <hydra/common.h>

namespace hydra {

inline double time_evolve_inplace(BondList const &bonds, StateCplx &state, double time,
                           double precision = 1e-12) {
  // auto const &block = state.block();

  // int iter = 1;
  // arma::cx_vec apply_A = [&iter, &bonds, &block](arma::cx_vec const &v) {
  //   auto ta = rightnow();
  //   apply(bonds, block, v, block, w);
  //   Log(2, "Lanczos iteration {}", iter);
  //   timing(ta, rightnow(), "MVM", 2);
  //   ++iter;
  // };
  // auto &v0 = state.vector();
  // auto t0 = rightnow();
  // auto [err, hump] = zahexpv(time, apply_A, v0, )
  // timing(t0, rightnow(), "Lanczos time", 1);
  return exp_sym_v_inplace(bonds, state, complex(0., -time), precision);
}

template <class StateT>
inline StateT time_evolve(BondList const &bonds, StateT state, double tau,
                   double precision = 1e-12) {
  time_evolve_inplace(bonds, state, tau, precision);
  return state;
}

template <class StateT>
inline double imag_time_evolve_inplace(BondList const &bonds, StateT &state,
                                double tau, double precision = 1e-12) {
  return exp_sym_v_inplace(bonds, state, -tau, precision);
}

template <class StateT>
inline StateT imag_time_evolve(BondList const &bonds, StateT state, double tau,
                        double precision = 1e-12) {
  imag_time_evolve_inplace(bonds, state, tau, precision);
  return state;
}

} // namespace hydra
