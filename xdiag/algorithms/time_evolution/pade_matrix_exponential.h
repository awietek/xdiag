//
// Created by Luke Staszewski on 07.02.23.
//
#pragma once

#include <algorithm>
#include <cmath>
#include <xdiag/extern/armadillo/armadillo>

namespace xdiag {

template <typename coeff_t>
arma::Mat<coeff_t> expm(arma::Mat<coeff_t> const &A, coeff_t alpha = 1.);

} // namespace xdiag
