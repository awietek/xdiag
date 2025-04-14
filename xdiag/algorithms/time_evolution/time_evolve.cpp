// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "time_evolve.hpp"

#include <xdiag/algorithms/time_evolution/evolve_lanczos.hpp>
#include <xdiag/algorithms/time_evolution/time_evolve_expokit.hpp>

namespace xdiag {

State time_evolve(OpSum const &H, State psi, double time, double precision,
                  std::string algorithm) try {
  time_evolve_inplace(H, psi, time, precision, algorithm);
  return psi;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void time_evolve_inplace(OpSum const &H, State &psi, double time,
                         double precision, std::string algorithm) try {
  if (algorithm == "lanczos") {
    // minus sign in exp(-iHt) implemented here
    evolve_lanczos_inplace(H, psi, complex(0, -time), precision);
  } else if (algorithm == "expokit") {
    // minus sign in exp(-iHt) implemented in expokit routine
    time_evolve_expokit_inplace(H, psi, time, precision);
  } else {
    XDIAG_THROW(
        fmt::format("Invalid time-evolution algorithm specified: \"{}\". Must "
                    "be one of \"lanczos\" or \"expokit\"."));
  }

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
