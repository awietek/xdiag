// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

#include <xdiag/extern/armadillo/armadillo>

namespace xdiag {

template <typename coeff_t>
arma::Mat<coeff_t> read_vectors(std::string type, std::string path_to_vecs,
                                int n = 0);
arma::mat read_vectors_real(std::string type, std::string path_to_vecs,
                            int n = 0);
arma::cx_mat read_vectors_cplx(std::string type, std::string path_to_vecs,
                               int n = 0);

} // namespace xdiag
