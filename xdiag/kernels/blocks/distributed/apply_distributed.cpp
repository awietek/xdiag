// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply_distributed.hpp"

#include <xdiag/utils/variants.hpp>

namespace xdiag::kernels {

template <typename coeff_t>
void apply_distributed(OpSum const &ops, Block const &block_in,
                       arma::Mat<coeff_t> const &mat_in, Block const &block_out,
                       arma::Mat<coeff_t> &mat_out) {

  utils::visit_same_type(
      block_in, block_out,
      [&](auto &&bin, auto &&bout) {
        if (mat_in.n_cols != mat_out.n_cols) {
          XDIAG_THROW("Input and output matrix do not have the same number "
                      "of columns");
        }

        // Apply column-wise
        int64_t ncols = mat_in.n_cols;
        for (int i = 0; i < ncols; ++i) {
          arma::Col<coeff_t> vec_in(const_cast<coeff_t *>(mat_in.colptr(i)),
                                    mat_in.n_rows, false, true);
          arma::Col<coeff_t> vec_out(mat_out.colptr(i), mat_out.n_rows, false,
                                     true);
          apply_distributed(ops, basis_in, vec_in, basis_out, vec_out);
        }
      },
      "Input and output block are not of the same type");
}

template void apply_distributed(OpSum const &, Block const &block_in,
                                arma::mat const &, Block const &block_out,
                                arma::mat &);
template void apply_distributed(OpSum const &, Block const &block_in,
                                arma::cx_mat const &, Block const &block_out,
                                arma::cx_mat &);
} // namespace xdiag::kernels
