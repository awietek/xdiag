// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

// Builds the CSR-style arrays (ptr of length ndim+1, idx, data) for a block.
//   transpose == false -> CSR: ndim = nrows, groups are rows,    idx = columns
//   transpose == true  -> CSC: ndim = ncols, groups are columns, idx = rows
// This is the single place the (expensive) per-basis-type dispatch is
// instantiated; csr_matrix.cpp and csc_matrix.cpp both call it so the visitor
// over all concrete basis types is compiled only once rather than in each.
template <typename idx_t, typename coeff_t, typename block_t>
void build_csr_arrays(OpSum const &ops, block_t const &block_in,
                      block_t const &block_out, idx_t ndim, idx_t i0,
                      bool transpose, arma::Col<idx_t> &ptr,
                      arma::Col<idx_t> &idx, arma::Col<coeff_t> &data);

} // namespace xdiag
