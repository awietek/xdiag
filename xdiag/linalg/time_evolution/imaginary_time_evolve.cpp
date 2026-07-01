// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "imaginary_time_evolve.hpp"

#include <xdiag/linalg/time_evolution/evolve_lanczos.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

template <typename op_t>
static State imaginary_time_evolve(op_t const &H, State psi, double time,
                                   double precision, double shift) try {
  imaginary_time_evolve_inplace(H, psi, time, precision, shift);
  return psi;
}
XDIAG_CATCH

State imaginary_time_evolve(OpSum const &H, State psi, double time,
                            double precision, double shift) try {
  return imaginary_time_evolve<OpSum>(H, psi, time, precision, shift);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
State imaginary_time_evolve(CSRMatrix<idx_t, coeff_t> const &H, State psi,
                            double time, double precision, double shift) try {
  return imaginary_time_evolve<CSRMatrix<idx_t, coeff_t>>(H, psi, time,
                                                          precision, shift);
}
XDIAG_CATCH

#define XDIAG_INST(IDX, COEFF)                                                 \
  template State imaginary_time_evolve(CSRMatrix<IDX, COEFF> const &, State,   \
                                       double, double, double);
XDIAG_INST(int32_t, double)
XDIAG_INST(int32_t, complex)
XDIAG_INST(int64_t, double)
XDIAG_INST(int64_t, complex)
#undef XDIAG_INST

template <typename op_t>
static void imaginary_time_evolve_inplace(op_t const &H, State &psi,
                                          double time, double precision,
                                          double shift) try {
  // minus sign in exp(-Ht) implemented here
  evolve_lanczos_inplace(H, psi, -time, precision, shift);
}
XDIAG_CATCH

void imaginary_time_evolve_inplace(OpSum const &H, State &psi, double time,
                                   double precision, double shift) try {
  imaginary_time_evolve_inplace<OpSum>(H, psi, time, precision, shift);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
void imaginary_time_evolve_inplace(CSRMatrix<idx_t, coeff_t> const &H,
                                   State &psi, double time, double precision,
                                   double shift) try {
  imaginary_time_evolve_inplace<CSRMatrix<idx_t, coeff_t>>(H, psi, time,
                                                           precision, shift);
}
XDIAG_CATCH

#define XDIAG_INST(IDX, COEFF)                                                 \
  template void imaginary_time_evolve_inplace(CSRMatrix<IDX, COEFF> const &,   \
                                              State &, double, double, double);
XDIAG_INST(int32_t, double)
XDIAG_INST(int32_t, complex)
XDIAG_INST(int64_t, double)
XDIAG_INST(int64_t, complex)
#undef XDIAG_INST

} // namespace xdiag
