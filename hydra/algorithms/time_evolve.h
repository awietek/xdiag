#pragma once

#include <hydra/algorithms/exp_sym_v.h>
#include <hydra/common.h>

namespace hydra {

template <class StateT>
double time_evolve_inplace(BondList const &bonds, StateT &state, double tau,
                           double precision = 1e-12) {
  return exp_sym_v_inplace(bonds, state, complex(0., -tau), precision);
}

template <class StateT>
StateT time_evolve(BondList const &bonds, StateT state, double tau,
                   double precision = 1e-12) {
  time_evolve_inplace(bonds, state, tau, precision);
  return state;
}

template <class StateT>
double imag_time_evolve_inplace(BondList const &bonds, StateT &state,
                                double tau, double precision = 1e-12) {
  return exp_sym_v_inplace(bonds, state, -tau, precision);
}

template <class StateT>
StateT imag_time_evolve(BondList const &bonds, StateT state, double tau,
                        double precision = 1e-12) {
  imag_time_evolve_inplace(bonds, state, tau, precision);
  return state;
}

} // namespace hydra
