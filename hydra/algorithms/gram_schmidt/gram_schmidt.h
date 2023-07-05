#pragma once

#include <string>

#include <extern/armadillo/armadillo>
#include <hydra/common.h>

namespace hydra {

template <typename coeff_t>
void gram_schmidt_inplace(arma::Mat<coeff_t> &M, int64_t max_col = -1,
                          std::string method = "mgs", int iterations = 2);

template <typename coeff_t>
arma::Mat<coeff_t> gram_schmidt(arma::Mat<coeff_t> const &M,
                                int64_t max_col = -1,
                                std::string method = "mgs", int iterations = 2);

} // namespace hydra
