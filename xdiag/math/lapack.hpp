// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <limits>
#include <type_traits>
#include <vector>

#include <xdiag/armadillo.hpp>
#include <xdiag/utils/error.hpp>

// LAPACK routines not wrapped by armadillo. Symbol naming and integer width
// follow the armadillo configuration set in xdiag/armadillo.hpp.

#if defined(ARMA_BLAS_CAPITALS)
#define XDIAG_LAPACK_ssterf SSTERF
#define XDIAG_LAPACK_dsterf DSTERF
#define XDIAG_LAPACK_sstebz SSTEBZ
#define XDIAG_LAPACK_dstebz DSTEBZ
#else
#define XDIAG_LAPACK_ssterf ssterf
#define XDIAG_LAPACK_dsterf dsterf
#define XDIAG_LAPACK_sstebz sstebz
#define XDIAG_LAPACK_dstebz dstebz
#endif

#if defined(ARMA_BLAS_UNDERSCORE)
#define XDIAG_FORTRAN_B(function) function##_
#else
#define XDIAG_FORTRAN_B(function) function
#endif
#define XDIAG_FORTRAN(function) XDIAG_FORTRAN_B(function)

// Trailing lengths fortran passes for each character argument.
#if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
#define XDIAG_FORTRAN_CHARLEN2_DECL , arma::blas_len len1, arma::blas_len len2
#define XDIAG_FORTRAN_CHARLEN2_ARGS , 1, 1
#else
#define XDIAG_FORTRAN_CHARLEN2_DECL
#define XDIAG_FORTRAN_CHARLEN2_ARGS
#endif

extern "C" {
void XDIAG_FORTRAN(XDIAG_LAPACK_ssterf)(const arma::blas_int *n, float *d,
                                        float *e, arma::blas_int *info);
void XDIAG_FORTRAN(XDIAG_LAPACK_dsterf)(const arma::blas_int *n, double *d,
                                        double *e, arma::blas_int *info);
void XDIAG_FORTRAN(XDIAG_LAPACK_sstebz)(
    const char *range, const char *order, const arma::blas_int *n,
    const float *vl, const float *vu, const arma::blas_int *il,
    const arma::blas_int *iu, const float *abstol, const float *d,
    const float *e, arma::blas_int *m, arma::blas_int *nsplit, float *w,
    arma::blas_int *iblock, arma::blas_int *isplit, float *work,
    arma::blas_int *iwork, arma::blas_int *info XDIAG_FORTRAN_CHARLEN2_DECL);
void XDIAG_FORTRAN(XDIAG_LAPACK_dstebz)(
    const char *range, const char *order, const arma::blas_int *n,
    const double *vl, const double *vu, const arma::blas_int *il,
    const arma::blas_int *iu, const double *abstol, const double *d,
    const double *e, arma::blas_int *m, arma::blas_int *nsplit, double *w,
    arma::blas_int *iblock, arma::blas_int *isplit, double *work,
    arma::blas_int *iwork, arma::blas_int *info XDIAG_FORTRAN_CHARLEN2_DECL);
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

// Eigenvalues il to iu (one based, ascending) of a symmetric tridiagonal
// matrix by bisection, written to eigvals. The n diagonal and n-1 offdiagonal
// entries are left unchanged.
template <typename T>
inline void stebz(int64_t n, T const *diag, T const *offdiag, int64_t il,
                  int64_t iu, T *eigvals) try {
  static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                "stebz is only defined for float and double");
  if (n > std::numeric_limits<arma::blas_int>::max()) {
    XDIAG_THROW("dimension too large for the LAPACK integer type");
  }
  if ((il < 1) || (iu < il) || (iu > n)) {
    XDIAG_THROW("invalid eigenvalue index range");
  }
  arma::blas_int nn = n, iil = il, iiu = iu, m = 0, nsplit = 0, info = 0;
  T vl = 0, vu = 0, abstol = 0;
  std::vector<T> w(n), work(4 * n);
  std::vector<arma::blas_int> iblock(n), isplit(n), iwork(3 * n);
  if constexpr (std::is_same_v<T, float>) {
    XDIAG_FORTRAN(XDIAG_LAPACK_sstebz)(
        "I", "E", &nn, &vl, &vu, &iil, &iiu, &abstol, diag, offdiag, &m,
        &nsplit, w.data(), iblock.data(), isplit.data(), work.data(),
        iwork.data(), &info XDIAG_FORTRAN_CHARLEN2_ARGS);
  } else {
    XDIAG_FORTRAN(XDIAG_LAPACK_dstebz)(
        "I", "E", &nn, &vl, &vu, &iil, &iiu, &abstol, diag, offdiag, &m,
        &nsplit, w.data(), iblock.data(), isplit.data(), work.data(),
        iwork.data(), &info XDIAG_FORTRAN_CHARLEN2_ARGS);
  }
  if ((info != 0) || (m != iiu - iil + 1)) {
    XDIAG_THROW("LAPACK stebz did not return the requested eigenvalues");
  }
  std::copy(w.begin(), w.begin() + m, eigvals);
}
XDIAG_CATCH

} // namespace xdiag::math::lapack

#undef XDIAG_LAPACK_ssterf
#undef XDIAG_LAPACK_dsterf
#undef XDIAG_LAPACK_sstebz
#undef XDIAG_LAPACK_dstebz
#undef XDIAG_FORTRAN_B
#undef XDIAG_FORTRAN
#undef XDIAG_FORTRAN_CHARLEN2_DECL
#undef XDIAG_FORTRAN_CHARLEN2_ARGS
