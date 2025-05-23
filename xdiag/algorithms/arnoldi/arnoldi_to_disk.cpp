// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "arnoldi_to_disk.hpp"

#include <xdiag/extern/fmt/format.hpp>
#include <xdiag/common.hpp>
#include <xdiag/utils/read_vectors.hpp>

namespace xdiag {

template <typename coeff_t>
arma::Mat<coeff_t> read_arnoldi_vectors(std::string path_to_vecs, int n) {
  return read_vectors<coeff_t>("Arnoldi", path_to_vecs, n);
}

template arma::Mat<double>
read_arnoldi_vectors<double>(std::string path_to_vecs, int n);
template arma::Mat<complex>
read_arnoldi_vectors<complex>(std::string path_to_vecs, int n);

arma::mat read_arnoldi_vectors_real(std::string path_to_vecs, int n) {
  return read_arnoldi_vectors<double>(path_to_vecs, n);
}
arma::cx_mat read_arnoldi_vectors_cplx(std::string path_to_vecs, int n) {
  return read_arnoldi_vectors<complex>(path_to_vecs, n);
}

template <typename coeff_t>
arma::Mat<coeff_t> read_ritz_vectors(std::string path_to_vecs, int n) {
  return read_vectors<coeff_t>("Ritz", path_to_vecs, n);
}

template arma::Mat<double> read_ritz_vectors<double>(std::string path_to_vecs,
                                                     int n);
template arma::Mat<complex> read_ritz_vectors<complex>(std::string path_to_vecs,
                                                       int n);

arma::mat read_ritz_vectors_real(std::string path_to_vecs, int n) {
  return read_ritz_vectors<double>(path_to_vecs, n);
}
arma::cx_mat read_ritz_vectors_cplx(std::string path_to_vecs, int n) {
  return read_ritz_vectors<complex>(path_to_vecs, n);
}

} // namespace xdiag
