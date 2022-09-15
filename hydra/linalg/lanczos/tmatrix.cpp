#include "tmatrix.h"

namespace hydra {

void Tmatrix::append(double alpha, double beta) {
  alphas_.push_back(alpha);
  betas_.push_back(beta);
}

void Tmatrix::pop() {
  alphas_.pop_back();
  betas_.pop_back();
}

arma::vec Tmatrix::eigenvalues() const {
  if (size() == 0) {
    return arma::Col<double>();
  } else {
    int N = alphas_.size();
    assert(N == (int)betas_.size());
    auto bts = betas()(arma::span(0, N - 1));
    auto tmat = arma::diagmat(alphas()) + arma::diagmat(bts, 1) +
                arma::diagmat(bts, -1);

    arma::vec eigs;
    arma::eig_sym(eigs, tmat);
    return eigs;
  }
}

arma::mat Tmatrix::eigenvectors() const {
  if (size() == 0) {
    return arma::Mat<double>();
  } else {
    int N = alphas_.size();
    assert(N == (int)betas_.size());
    auto bts = betas()(arma::span(0, N - 1));
    auto tmat = arma::diagmat(alphas()) + arma::diagmat(bts, 1) +
                arma::diagmat(bts, -1);

    arma::vec eigs;
    arma::mat evecs;
    arma::eig_sym(eigs, evecs, tmat);
    return evecs;
  }
}

std::pair<arma::vec, arma::mat> Tmatrix::eigen() const {
  if (size() == 0) {
    return {arma::vec(), arma::mat()};
  } else {
    int N = alphas_.size();
    assert(N == (int)betas_.size());
    auto bts = betas()(arma::span(0, N - 1));
    auto tmat = arma::diagmat(alphas()) + arma::diagmat(bts, 1) +
                arma::diagmat(bts, -1);

    arma::vec eigs;
    arma::mat evecs;
    arma::eig_sym(eigs, evecs, tmat);
    return {eigs, evecs};
  }
}

} // namespace hydra
