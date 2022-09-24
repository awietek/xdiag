#pragma once

#include <utility>

#include "extern/armadillo/armadillo"

#include <hydra/common.h>

namespace hydra {

class Tmatrix {
public:
  Tmatrix() = default;
  Tmatrix(std::vector<double> const &alphas, std::vector<double> const &betas)
      : alphas_(alphas), betas_(betas) {
    assert(alphas.size() == betas.size());
  }
  void append(double alpha, double beta);
  void pop();

  idx_t size() const { return alphas_.size(); }

  arma::vec alphas() const { return arma::vec(alphas_); }
  arma::vec betas() const { return betas_; }
  arma::mat mat() const;

  arma::vec eigenvalues() const;
  arma::mat eigenvectors() const;
  std::pair<arma::vec, arma::mat> eigen() const;

  void print_log() const;

private:
  std::vector<double> alphas_;
  std::vector<double> betas_;
};

} // namespace hydra
