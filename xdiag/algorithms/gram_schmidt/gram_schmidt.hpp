#pragma once

#include <string>

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/common.hpp>

namespace xdiag {

template <typename coeff_t>
void gram_schmidt_inplace(arma::Mat<coeff_t> &M, int64_t max_col = -1,
                          int iterations = 2);

template <typename coeff_t>
arma::Mat<coeff_t> gram_schmidt(arma::Mat<coeff_t> const &M,
                                int64_t max_col = -1, int iterations = 2);

} // namespace xdiag
