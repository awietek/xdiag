#pragma once

#include <lila/all.h>
#include <utility>

namespace hydra {

template <class coeff_t> class Tmatrix {
public:
  using size_type = lila::size_type;
  using real_type = lila::real_t<coeff_t>;

  Tmatrix() = default;
  Tmatrix(lila::Vector<real_type> const &alphas,
                 lila::Vector<real_type> const &betas)
      : alphas_(alphas), betas_(betas) {
    assert(alphas.size() == betas.size());
  }
  void append(real_type alpha, real_type beta);

  size_type size() const {return alphas_.size(); }

  lila::Vector<real_type> const &alphas() const { return alphas_; }
  lila::Vector<real_type> const &betas() const { return betas_; }

  lila::Vector<real_type> eigenvalues() const;
  std::pair<lila::Vector<real_type>, lila::Matrix<real_type>> eigen() const;

private:
  lila::Vector<real_type> alphas_;
  lila::Vector<real_type> betas_;
};

} // namespace hydra
