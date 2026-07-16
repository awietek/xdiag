// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>

namespace xdiag::testcases {

// Robust hermiticity / symmetry check for test code. arma's Mat::is_hermitian
// is unreliable for small complex matrices: a ~1e-15 imaginary residue on the
// diagonal (e.g. from complex symmetry characters) makes it report false even
// though the matrix is hermitian to machine precision. We instead test the
// actual maximum deviation |A - A^dagger|. Works for real (symmetry) and complex
// (hermiticity) matrices alike.
template <typename T>
inline bool is_approx_hermitian(arma::Mat<T> const &A, double tol = 1e-8) {
  return arma::abs(A - A.t()).max() < tol;
}

} // namespace xdiag::testcases
