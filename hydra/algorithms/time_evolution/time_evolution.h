#pragma once

#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/states/state.h>

namespace hydra {

std::tuple<double, double> time_evolve_inplace(BondList const &bonds,
                                               State &state, double time,
                                               double precision = 1e-12,
                                               int64_t m = 5, double anorm = 0.,
                                               int64_t nnorm = 2);

State time_evolve(BondList const &bonds, State state, double time,
                  double precision = 1e-12, int64_t m = 5, double anorm = 0.,
                  int64_t nnorm = 2);

double imag_time_evolve_inplace(BondList const &bonds, State &state,
                                double time, double precision = 1e-12);

State imag_time_evolve(BondList const &bonds, State state, double time,
                       double precision = 1e-12);

} // namespace hydra
