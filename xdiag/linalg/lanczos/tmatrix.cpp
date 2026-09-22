// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "tmatrix.hpp"

#include <algorithm>
#include <cassert>
#include <string>

#include <xdiag/math/lapack.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag {

Tmatrix::Tmatrix(std::vector<double> const &alphas,
                 std::vector<double> const &betas)
    : alphas_(alphas), betas_(betas) {
  assert(alphas.size() == betas.size());
}

void Tmatrix::append(double alpha, double beta) {
  alphas_.push_back(alpha);
  betas_.push_back(beta);
}

void Tmatrix::pop() {
  alphas_.pop_back();
  betas_.pop_back();
}

arma::mat Tmatrix::mat() const try {
  if (size() == 0) {
    return arma::Mat<double>();
  } else if (size() == 1) {
    return arma::Mat<double>(1, 1, arma::fill::value(alphas_[0]));
  } else {
    int N = alphas_.size();
    assert(N == (int)betas_.size());
    arma::vec bts = betas().subvec(0, N - 2);
    arma::mat tmat = arma::diagmat(alphas()) + arma::diagmat(bts, 1) +
                     arma::diagmat(bts, -1);
    return tmat;
  }
} catch (...) {
  XDIAG_THROW("cannot armadillo matrix from Tmatrix");
  return arma::mat();
}

arma::vec Tmatrix::eigenvalues() const try {
  if (size() == 0) {
    return arma::Col<double>();
  } else if (size() == 1) {
    return arma::Col<double>(1, arma::fill::value(alphas_[0]));
  } else {
    int64_t n = size();
    // sterf overwrites both arrays, so the stored coefficients are copied
    arma::vec diag(alphas_.data(), n);
    arma::vec offdiag(betas_.data(), n - 1);
    math::lapack::sterf(n, diag.memptr(), offdiag.memptr());
    return diag;
  }
} catch (...) {
  XDIAG_THROW("cannot compute eigenvalues of Tmatrix");
  return arma::vec();
}

arma::vec Tmatrix::eigenvalues_lowest(int64_t k) const try {
  int64_t n = size();
  if ((n == 0) || (k >= n)) {
    return eigenvalues();
  }
  arma::vec eigs(k);
  math::lapack::stebz(n, alphas_.data(), betas_.data(), (int64_t)1, k,
                      eigs.memptr());
  return eigs;
} catch (...) {
  XDIAG_THROW("cannot compute lowest eigenvalues of Tmatrix");
  return arma::vec();
}

arma::mat Tmatrix::eigenvectors() const try {
  return eigen().second;
} catch (...) {
  XDIAG_THROW("cannot compute eigenvectors of Tmatrix");
  return arma::mat();
}

std::pair<arma::vec, arma::mat> Tmatrix::eigen() const try {
  if (size() == 0) {
    return {arma::vec(), arma::mat()};
  } else if (size() == 1) {
    return {arma::Col<double>(1, arma::fill::value(alphas_[0])),
            arma::Mat<double>(1, 1, arma::fill::value(1.0))};
  } else {
    int64_t n = size();
    // stedc overwrites diag with the eigenvalues and destroys offdiag
    arma::vec diag(alphas_.data(), n);
    arma::vec offdiag(betas_.data(), n - 1);
    arma::mat evecs(n, n);
    math::lapack::stedc(n, diag.memptr(), offdiag.memptr(), evecs.memptr());
    return {diag, evecs};
  }
} catch (...) {
  XDIAG_THROW("cannot compute eigendecomposition of Tmatrix");
  return std::pair<arma::vec, arma::mat>();
}

void Tmatrix::print_log() const try {
  arma::vec eigs = eigenvalues();
  double alpha = alphas_[size() - 1];
  double beta = betas_[size() - 1];
  Log(2, "alpha: {:.16f}", alpha);
  Log(2, "beta: {:.16f}", beta);
  std::string eigs_str;
  for (int64_t i = 0; i < std::min<int64_t>(3, eigs.size()); ++i) {
    eigs_str += fmt::format("{:.16f} ", eigs(i));
  }
  Log(2, "eigs: {}", eigs_str);
} catch (...) {
  XDIAG_THROW("cannot logging of Tmatrix");
}

bool Tmatrix::operator==(Tmatrix const &rhs) const {
  return (alphas_ == rhs.alphas_) && (betas_ == rhs.betas_);
}
bool Tmatrix::operator!=(Tmatrix const &rhs) const { return !operator==(rhs); }

} // namespace xdiag
