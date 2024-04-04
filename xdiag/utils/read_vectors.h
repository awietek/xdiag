#pragma once

#include <xdiag/extern/armadillo/armadillo>
#include <string>

namespace xdiag {

template <typename coeff_t>
arma::Mat<coeff_t> read_vectors(std::string type, std::string path_to_vecs,
                                int n = 0);
arma::mat read_vectors_real(std::string type, std::string path_to_vecs,
                            int n = 0);
arma::cx_mat read_vectors_cplx(std::string type, std::string path_to_vecs,
                               int n = 0);

} // namespace xdiag
