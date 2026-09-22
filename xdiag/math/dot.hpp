// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/math/complex.hpp>

namespace xdiag::math {

// Internal routines (used to differentiate dot product for distributed blocks)
double dot(Block const &block, arma::vec const &v, arma::vec const &w);
complex dot(Block const &block, arma::cx_vec const &v, arma::cx_vec const &w);

// Internal. Shared-memory dot product conj(v) . w whose result is bitwise
// independent of the number of OpenMP threads: the vectors are summed in
// fixed-length chunks and the chunk sums are reduced in order (see dot.cpp).
template <typename coeff_t>
coeff_t dot_chunked(arma::Col<coeff_t> const &v, arma::Col<coeff_t> const &w);

template <typename coeff_t>
arma::Mat<coeff_t> matrix_dot(Block const &block, arma::Mat<coeff_t> const &V,
                              arma::Mat<coeff_t> const &W);

} // namespace xdiag::math
