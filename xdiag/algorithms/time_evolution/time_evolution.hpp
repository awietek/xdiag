#pragma once

#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

std::tuple<double, double> time_evolve_inplace(OpSum const &ops, State &state,
                                               double time,
                                               double precision = 1e-12,
                                               int64_t m = 5, double anorm = 0.,
                                               int64_t nnorm = 2);

State time_evolve(OpSum const &ops, State state, double time,
                  double precision = 1e-12, int64_t m = 5, double anorm = 0.,
                  int64_t nnorm = 2);

double imag_time_evolve_inplace(OpSum const &ops, State &state, double time,
                                double precision = 1e-12);

State imag_time_evolve(OpSum const &ops, State state, double time,
                       double precision = 1e-12);

} // namespace xdiag
