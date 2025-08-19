// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <xdiag/algebra/sparse/sparse_matrix_types.hpp>
#include <xdiag/common.hpp>

namespace xdiag {

template <typename idx_t, typename coeff_t>
constexpr bool isreal(CSRMatrix<idx_t, coeff_t> const &) {
  return isreal<coeff_t>();
}
template <typename idx_t, typename coeff_t>
constexpr bool isreal(CSCMatrix<idx_t, coeff_t> const &) {
  return isreal<coeff_t>();
}
template <typename idx_t, typename coeff_t>
constexpr bool isreal(COOMatrix<idx_t, coeff_t> const &) {
  return isreal<coeff_t>();
}

template <typename idx_t, typename coeff_t>
inline bool ishermitian(CSRMatrix<idx_t, coeff_t> const &A) {
  return A.ishermitian;
}
template <typename idx_t, typename coeff_t>
inline bool ishermitian(CSCMatrix<idx_t, coeff_t> const &A) {
  return A.ishermitian;
}
template <typename idx_t, typename coeff_t>
inline bool ishermitian(COOMatrix<idx_t, coeff_t> const &A) {
  return A.ishermitian;
}

} // namespace xdiag
