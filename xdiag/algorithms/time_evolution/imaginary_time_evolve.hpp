// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>

#include <xdiag/algebra/sparse/sparse_matrix_types.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

XDIAG_API State imaginary_time_evolve(OpSum const &H, State psi, double time,
                                      double precision = 1e-12,
                                      double shift = 0.);
template <typename idx_t, typename coeff_t>
XDIAG_API State imaginary_time_evolve(CSRMatrix<idx_t, coeff_t> const &H,
                                      State psi, double time,
                                      double precision = 1e-12);

XDIAG_API void imaginary_time_evolve_inplace(OpSum const &H, State &psi,
                                             double time,
                                             double precision = 1e-12,
                                             double shift = 0.);

template <typename idx_t, typename coeff_t>
XDIAG_API void imaginary_time_evolve_inplace(CSRMatrix<idx_t, coeff_t> const &H,
                                             State &psi, double time,
                                             double precision = 1e-12);

} // namespace xdiag
