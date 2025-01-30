#pragma once

#include <xdiag/common.hpp>

#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

XDIAG_API State imaginary_time_evolve(OpSum const &H, State psi, double time,
                                      double precision = 1e-12,
                                      double shift = 0.);

XDIAG_API void imaginary_time_evolve_inplace(OpSum const &H, State &psi,
                                             double time,
                                             double precision = 1e-12,
					     double shift = 0.);

} // namespace xdiag
