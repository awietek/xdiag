// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

// We stick to armadillo version 12.8.4 to minimize compiler constraints

#define ARMA_DONT_USE_WRAPPER

#ifdef XDIAG_USE_MKL
#define ARMA_64BIT_WORD
#define ARMA_BLAS_LONG_LONG
#else
#define ARMA_64BIT_WORD
#define ARMA_BLAS_UNDERSCORE
#endif

#ifdef XDIAG_USE_HDF5
// #define ARMA_USE_HDF5
// commented out since caused an error
// error: no matching function for call to 'H5Ovisit3'
// on OSX
#endif

// Redefinition of "Op" to avoid name colission with xdiag::Op
#define Op XDIAG_OP_SUBSTITUTE
#include <extern/armadillo/armadillo>
#undef Op
