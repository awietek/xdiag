#pragma once

#include <utility>

#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/common.hpp>

namespace xdiag {

class Tmatrix {
public:
  Tmatrix() = default;
  Tmatrix(std::vector<double> const &alphas, std::vector<double> const &betas);
  void append(double alpha, double beta);
  void pop();

  int64_t size() const { return alphas_.size(); }

  arma::vec alphas() const { return arma::vec(alphas_); }
  arma::vec betas() const { return betas_; }
  arma::mat mat() const;

  arma::vec eigenvalues() const;
  arma::mat eigenvectors() const;
  std::pair<arma::vec, arma::mat> eigen() const;

  void print_log() const;

  bool operator==(Tmatrix const& rhs) const;
  bool operator!=(Tmatrix const& rhs) const;
  
private:
  std::vector<double> alphas_;
  std::vector<double> betas_;
};

} // namespace xdiag
