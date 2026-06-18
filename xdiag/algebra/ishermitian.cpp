// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "ishermitian.hpp"

#include <cstdint>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/operators/hc.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

bool ishermitian(OpSum const &ops, algebra::Algebra const &algebra,
                 double tol) try {
  return isapprox(ops, hc(ops), algebra, tol, tol);
}
XDIAG_CATCH

bool ishermitian(OpSum const &ops, Block const &block, double tol) try {
  return ishermitian(ops, algebra::symmetry_algebra(block), tol);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
bool ishermitian(CSRMatrix<idx_t, coeff_t> const &A, Block const &, double) {
  return A.ishermitian;
}
template <typename idx_t, typename coeff_t>
bool ishermitian(CSCMatrix<idx_t, coeff_t> const &A, Block const &, double) {
  return A.ishermitian;
}
template <typename idx_t, typename coeff_t>
bool ishermitian(COOMatrix<idx_t, coeff_t> const &A, Block const &, double) {
  return A.ishermitian;
}

// Explicit instantiations matching the sparse-matrix coefficient/index types.
#define XDIAG_INSTANTIATE_ISHERMITIAN(Mat, idx_t, coeff_t)                     \
  template bool ishermitian(Mat<idx_t, coeff_t> const &, Block const &, double);

#define XDIAG_INSTANTIATE_ISHERMITIAN_ALL(Mat)                                 \
  XDIAG_INSTANTIATE_ISHERMITIAN(Mat, int32_t, double)                          \
  XDIAG_INSTANTIATE_ISHERMITIAN(Mat, int32_t, complex)                         \
  XDIAG_INSTANTIATE_ISHERMITIAN(Mat, int64_t, double)                          \
  XDIAG_INSTANTIATE_ISHERMITIAN(Mat, int64_t, complex)

XDIAG_INSTANTIATE_ISHERMITIAN_ALL(CSRMatrix)
XDIAG_INSTANTIATE_ISHERMITIAN_ALL(CSCMatrix)
XDIAG_INSTANTIATE_ISHERMITIAN_ALL(COOMatrix)

#undef XDIAG_INSTANTIATE_ISHERMITIAN_ALL
#undef XDIAG_INSTANTIATE_ISHERMITIAN

} // namespace xdiag
