// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "arma_to_cx.hpp"

namespace xdiag::utils {
arma::cx_vec to_cx_vec(arma::vec const &A) {
  return arma::cx_vec(A, arma::vec(A.n_rows, A.n_cols, arma::fill::zeros));
}

arma::cx_vec to_cx_vec(arma::cx_vec const &A) { return A; }

arma::cx_mat to_cx_mat(arma::mat const &A) {
  return arma::cx_mat(A, arma::mat(A.n_rows, A.n_cols, arma::fill::zeros));
}

arma::cx_mat to_cx_mat(arma::cx_mat const &A) { return A; }

} // namespace xdiag::utils
