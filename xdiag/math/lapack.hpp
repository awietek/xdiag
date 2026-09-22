// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <limits>
#include <type_traits>

#include <xdiag/armadillo.hpp>
#include <xdiag/utils/error.hpp>

// LAPACK routines not wrapped by armadillo. Symbol naming and integer width
// follow the armadillo configuration set in xdiag/armadillo.hpp.

#if defined(ARMA_BLAS_CAPITALS)
#define XDIAG_LAPACK_ssterf SSTERF
#define XDIAG_LAPACK_dsterf DSTERF
#else
#define XDIAG_LAPACK_ssterf ssterf
#define XDIAG_LAPACK_dsterf dsterf
#endif

#if defined(ARMA_BLAS_UNDERSCORE)
#define XDIAG_FORTRAN_B(function) function##_
#else
#define XDIAG_FORTRAN_B(function) function
#endif
#define XDIAG_FORTRAN(function) XDIAG_FORTRAN_B(function)

extern "C" {
void XDIAG_FORTRAN(XDIAG_LAPACK_ssterf)(const arma::blas_int *n, float *d,
                                        float *e, arma::blas_int *info);
void XDIAG_FORTRAN(XDIAG_LAPACK_dsterf)(const arma::blas_int *n, double *d,
                                        double *e, arma::blas_int *info);
}

namespace xdiag::math::lapack {

// Eigenvalues of a symmetric tridiagonal matrix with n diagonal entries and
// n-1 offdiagonal entries. Both arrays are overwritten, diag with the
// eigenvalues in ascending order.
template <typename T>
inline void sterf(int64_t n, T *diag, T *offdiag) try {
  static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                "sterf is only defined for float and double");
  if (n > std::numeric_limits<arma::blas_int>::max()) {
    XDIAG_THROW("dimension too large for the LAPACK integer type");
  }
  arma::blas_int nn = n;
  arma::blas_int info = 0;
  if constexpr (std::is_same_v<T, float>) {
    XDIAG_FORTRAN(XDIAG_LAPACK_ssterf)(&nn, diag, offdiag, &info);
  } else {
    XDIAG_FORTRAN(XDIAG_LAPACK_dsterf)(&nn, diag, offdiag, &info);
  }
  if (info != 0) {
    XDIAG_THROW("LAPACK sterf did not converge");
  }
}
XDIAG_CATCH

} // namespace xdiag::math::lapack

#undef XDIAG_LAPACK_ssterf
#undef XDIAG_LAPACK_dsterf
#undef XDIAG_FORTRAN_B
#undef XDIAG_FORTRAN
