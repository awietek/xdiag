#pragma once

#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

struct time_evolve_expokit_return_t {
  double error;
  double hump;
  State state;
};

XDIAG_API time_evolve_expokit_return_t time_evolve_expokit(
    OpSum const &H, State psi0, double time, double precision = 1e-12,
    int64_t m = 30, double anorm = 0., int64_t nnorm = 2);

struct time_evolve_expokit_inplace_return_t {
  double error;
  double hump;
};

XDIAG_API time_evolve_expokit_inplace_return_t time_evolve_expokit_inplace(
    OpSum const &H, State &psi, double time, double precision = 1e-12,
    int64_t m = 30, double anorm = 0., int64_t nnorm = 2);

} // namespace xdiag
