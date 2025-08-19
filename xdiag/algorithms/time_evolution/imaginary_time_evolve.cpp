// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "imaginary_time_evolve.hpp"

#include <xdiag/algorithms/time_evolution/evolve_lanczos.hpp>

namespace xdiag {

template <typename op_t>
static State imaginary_time_evolve(op_t const &H, State psi, double time,
                                   double precision, double shift) try {
  imaginary_time_evolve_inplace(H, psi, time, precision, shift);
  return psi;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

State imaginary_time_evolve(OpSum const &H, State psi, double time,
                            double precision, double shift) try {
  return imaginary_time_evolve<OpSum>(H, psi, time, precision, shift);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename idx_t, typename coeff_t>
State imaginary_time_evolve(CSRMatrix<idx_t, coeff_t> const &H, State psi,
                            double time, double precision, double shift) try {
  return imaginary_time_evolve<CSRMatrix<idx_t, coeff_t>>(H, psi, time,
                                                          precision, shift);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template State imaginary_time_evolve(CSRMatrix<int32_t, double> const &, State,
                                     double, double, double);
template State imaginary_time_evolve(CSRMatrix<int32_t, complex> const &, State,
                                     double, double, double);
template State imaginary_time_evolve(CSRMatrix<int64_t, double> const &, State,
                                     double, double, double);
template State imaginary_time_evolve(CSRMatrix<int64_t, complex> const &, State,
                                     double, double, double);

template <typename op_t>
static void imaginary_time_evolve_inplace(op_t const &H, State &psi,
                                          double time, double precision,
                                          double shift) try {
  // minus sign in exp(-Ht) implemented here
  evolve_lanczos_inplace(H, psi, -time, precision, shift);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void imaginary_time_evolve_inplace(OpSum const &H, State &psi, double time,
                                   double precision, double shift) try {
  imaginary_time_evolve_inplace<OpSum>(H, psi, time, precision, shift);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename idx_t, typename coeff_t>
void imaginary_time_evolve_inplace(CSRMatrix<idx_t, coeff_t> const &H,
                                   State &psi, double time, double precision,
                                   double shift) try {
  imaginary_time_evolve_inplace<CSRMatrix<idx_t, coeff_t>>(H, psi, time,
                                                           precision, shift);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void imaginary_time_evolve_inplace(CSRMatrix<int32_t, double> const &,
                                            State &, double, double, double);
template void imaginary_time_evolve_inplace(CSRMatrix<int32_t, complex> const &,
                                            State &, double, double, double);
template void imaginary_time_evolve_inplace(CSRMatrix<int64_t, double> const &,
                                            State &, double, double, double);
template void imaginary_time_evolve_inplace(CSRMatrix<int64_t, complex> const &,
                                            State &, double, double, double);

} // namespace xdiag
