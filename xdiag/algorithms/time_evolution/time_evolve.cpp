// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "time_evolve.hpp"

#include <xdiag/algorithms/time_evolution/evolve_lanczos.hpp>
#include <xdiag/algorithms/time_evolution/time_evolve_expokit.hpp>

namespace xdiag {

template <typename op_t>
static State time_evolve(op_t const &H, State psi, double time,
                         double precision, std::string algorithm) try {
  time_evolve_inplace(H, psi, time, precision, algorithm);
  return psi;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

State time_evolve(OpSum const &H, State psi, double time, double precision,
                  std::string algorithm) try {
  return time_evolve<OpSum>(H, psi, time, precision, algorithm);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename idx_t, typename coeff_t>
State time_evolve(CSRMatrix<idx_t, coeff_t> const &H, State psi, double time,
                  double precision, std::string algorithm) try {
  return time_evolve<CSRMatrix<idx_t, coeff_t>>(H, psi, time, precision,
                                                algorithm);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template State time_evolve(CSRMatrix<int32_t, double> const &, State, double,
                           double, std::string);
template State time_evolve(CSRMatrix<int32_t, complex> const &, State, double,
                           double, std::string);
template State time_evolve(CSRMatrix<int64_t, double> const &, State, double,
                           double, std::string);
template State time_evolve(CSRMatrix<int64_t, complex> const &, State, double,
                           double, std::string);

template <typename op_t>
static void time_evolve_inplace(op_t const &H, State &psi, double time,
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

void time_evolve_inplace(OpSum const &H, State &psi, double time,
                         double precision, std::string algorithm) try {
  time_evolve_inplace<OpSum>(H, psi, time, precision, algorithm);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename idx_t, typename coeff_t>
void time_evolve_inplace(CSRMatrix<idx_t, coeff_t> const &H, State &psi,
                         double time, double precision,
                         std::string algorithm) try {
  time_evolve_inplace<CSRMatrix<idx_t, coeff_t>>(H, psi, time, precision,
                                                 algorithm);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template void time_evolve_inplace(CSRMatrix<int32_t, double> const &, State &,
                                  double, double, std::string);
template void time_evolve_inplace(CSRMatrix<int32_t, complex> const &, State &,
                                  double, double, std::string);
template void time_evolve_inplace(CSRMatrix<int64_t, double> const &, State &,
                                  double, double, std::string);
template void time_evolve_inplace(CSRMatrix<int64_t, complex> const &, State &,
                                  double, double, std::string);
} // namespace xdiag
