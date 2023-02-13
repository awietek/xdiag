#pragma once

#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/states/state.h>

namespace hydra {

std::tuple<double, double> time_evolve_inplace(BondList const &bonds,
                                               StateCplx &state, double time,
                                               double precision = 1e-12,
                                               int m = 5, double anorm = 0.,
                                               int nnorm = 2);

template <typename coeff_t>
StateCplx time_evolve(BondList const &bonds, State<coeff_t> state, double time,
                      double precision = 1e-12, int m = 5, double anorm = 0.,
                      int nnorm = 2);

template <typename coeff_t>
double imag_time_evolve_inplace(BondList const &bonds, State<coeff_t> &state,
                                double time, double precision = 1e-12);

template <typename coeff_t>
State<coeff_t> imag_time_evolve(BondList const &bonds, State<coeff_t> state,
                                double time, double precision = 1e-12);

} // namespace hydra
