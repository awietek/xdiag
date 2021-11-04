#pragma once

#include <lila/all.h>
#include <utility>

namespace hydra {

class Tmatrix {
public:
  using size_type = lila::size_type;

  Tmatrix() = default;
  Tmatrix(lila::Vector<double> const &alphas,
	  lila::Vector<double> const &betas)
      : alphas_(alphas), betas_(betas) {
    assert(alphas.size() == betas.size());
  }
  void append(double alpha, double beta);
  void pop();

  size_type size() const {return alphas_.size(); }

  lila::Vector<double> const &alphas() const { return alphas_; }
  lila::Vector<double> const &betas() const { return betas_; }

  lila::Vector<double> eigenvalues() const;
  std::pair<lila::Vector<double>, lila::Matrix<double>> eigen() const;

private:
  lila::Vector<double> alphas_;
  lila::Vector<double> betas_;
};

} // namespace hydra
