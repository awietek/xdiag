#include "tmatrix.hpp"

#include <cassert>

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
    arma::vec eigs;
    arma::eig_sym(eigs, mat());
    return eigs;
  }
} catch (...) {
  XDIAG_THROW("cannot compute eigenvalues of Tmatrix");
  return arma::vec();
}

arma::mat Tmatrix::eigenvectors() const try {
  if (size() == 0) {
    return arma::Mat<double>();
  } else if (size() == 1) {
    return arma::Mat<double>(1, 1, arma::fill::value(1.0));
  } else {
    arma::vec eigs;
    arma::mat evecs;
    arma::eig_sym(eigs, evecs, mat());
    return evecs;
  }
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
    arma::vec eigs;
    arma::mat evecs;
    arma::eig_sym(eigs, evecs, mat());
    return {eigs, evecs};
  }
} catch (...) {
  XDIAG_THROW("cannot compute eigendecomposition of Tmatrix");
  return std::pair<arma::vec, arma::mat>();
}

void Tmatrix::print_log() const try {
  auto eigs = eigenvalues();
  double alpha = alphas_[size() - 1];
  double beta = betas_[size() - 1];
  Log(2, "alpha: {:.16f}", alpha);
  Log(2, "beta: {:.16f}", beta);
  if (eigs.size() == 1) {
    Log(2, "eigs: {:.16f}", eigs(0));
  } else if (eigs.size() == 2) {
    Log(2, "eigs: {:.16f} {:.16f}", eigs(0), eigs(1));
  } else {
    Log(2, "eigs: {:.16f} {:.16f} {:.16f}", eigs(0), eigs(1), eigs(2));
  }
} catch (...) {
  XDIAG_THROW("cannot logging of Tmatrix");
}

bool Tmatrix::operator==(Tmatrix const &rhs) const {
  return (alphas_ == rhs.alphas_) && (betas_ == rhs.betas_);
}
bool Tmatrix::operator!=(Tmatrix const &rhs) const { return !operator==(rhs); }

} // namespace xdiag
