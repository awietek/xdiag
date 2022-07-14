#pragma once

#include <utility>

#include <lila/all.h>

#include <hydra/common.h>

namespace hydra {

class Tmatrix {
public:
  Tmatrix() = default;
  Tmatrix(lila::Vector<double> const &alphas, lila::Vector<double> const &betas)
      : alphas_(alphas), betas_(betas) {
    assert(alphas.size() == betas.size());
  }
  void append(double alpha, double beta);
  void pop();

  idx_t size() const { return alphas_.size(); }

  lila::Vector<double> const &alphas() const { return alphas_; }
  lila::Vector<double> const &betas() const { return betas_; }

  lila::Vector<double> eigenvalues() const;
  lila::Matrix<double> eigenvectors() const;
  std::pair<lila::Vector<double>, lila::Matrix<double>> eigen() const;

private:
  lila::Vector<double> alphas_;
  lila::Vector<double> betas_;
};

} // namespace hydra
