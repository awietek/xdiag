#pragma once

#include <xdiag/common.hpp>

#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

XDIAG_API State time_evolve(OpSum const &H, State psi, double time,
                            double precision = 1e-12,
                            std::string algorithm = "lanczos");

XDIAG_API void time_evolve_inplace(OpSum const &H, State &psi, double time,
                                   double precision = 1e-12,
                                   std::string algorithm = "lanczos");

} // namespace xdiag
