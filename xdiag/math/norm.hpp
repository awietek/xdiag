// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>

namespace xdiag::math {

// Internal routines (used to differentiate dot product for distributed blocks)
template <typename coeff_t>
double norm(Block const &block, arma::Col<coeff_t> const &v);

template <typename coeff_t>
double norm1(Block const &block, arma::Col<coeff_t> const &v);

template <typename coeff_t>
double norminf(Block const &block, arma::Col<coeff_t> const &v);

} // namespace xdiag::math
