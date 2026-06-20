// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "dot.hpp"

#include <xdiag/utils/error.hpp>

#ifdef XDIAG_DISTRIBUTED
#include <xdiag/mpi/allreduce.hpp>
#include <xdiag/mpi/cdot_distributed.hpp>
#endif

namespace xdiag::math {

double dot(Block const &block, arma::vec const &v, arma::vec const &w) try {
#ifdef XDIAG_DISTRIBUTED
  if (isdistributed(block)) {
    return cdot_distributed(v, w);
  } else {
#else
  (void)block;
#endif
    return arma::dot(v, w);
#ifdef XDIAG_DISTRIBUTED
  }
#endif
}
XDIAG_CATCH

complex dot(Block const &block, arma::cx_vec const &v,
            arma::cx_vec const &w) try {
#ifdef XDIAG_DISTRIBUTED
  if (isdistributed(block)) {
    return cdot_distributed(v, w);
  } else {
#else
  (void)block;
#endif
    return arma::cdot(v, w);
#ifdef XDIAG_DISTRIBUTED
  }
#endif
}
XDIAG_CATCH

template <typename coeff_t>
arma::Mat<coeff_t> matrix_dot(Block const &block, arma::Mat<coeff_t> const &V,
                              arma::Mat<coeff_t> const &W) try {
  if (V.n_rows != V.n_rows) {
    XDIAG_THROW("Input matrices do not have the same number of rows");
  }
#ifdef XDIAG_DISTRIBUTED
  if (isdistributed(block)) {
    int64_t L = V.n_rows;
    int64_t m = V.n_cols;
    int64_t n = W.n_cols;
    arma::Mat<coeff_t> result(m, n, arma::fill::zeros);
    for (int64_t i = 0; i < m; ++i) {
      for (int64_t j = 0; j < n; ++j) {
        arma::Col<coeff_t> cv(const_cast<coeff_t *>(V.colptr(i)), L, false,
                              true);
        arma::Col<coeff_t> cw(const_cast<coeff_t *>(W.colptr(j)), L, false,
                              true);
        result(i, j) = cdot_distributed(cv, cw);
      }
    }
    return result;
  } else {
#else
  (void)block;
#endif
    return V.t() * W;
#ifdef XDIAG_DISTRIBUTED
  }
#endif
}
XDIAG_CATCH

template arma::mat matrix_dot(Block const &, arma::mat const &,
                              arma::mat const &);
template arma::cx_mat matrix_dot(Block const &, arma::cx_mat const &,
                                 arma::cx_mat const &);

} // namespace xdiag::math
