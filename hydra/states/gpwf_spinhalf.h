#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/common.h>

namespace hydra {

template <class coeff_t> class GPWFSpinhalf {
public:
  GPWFSpinhalf(int n_sites, arma::Mat<coeff_t> onebody_wfs, int n_up = -1);
  int n_sites() const { return n_sites_; }
  int n_up() const { return n_up_; }
  int n_dn() const { return n_dn_; }

  coeff_t coefficient(uint64_t state);

  bool operator==(GPWFSpinhalf<coeff_t> const &other);

  arma::Mat<coeff_t> onebody_wfs_up() const { return onebody_wfs_up_; }
  arma::Mat<coeff_t> onebody_wfs_dn() const { return onebody_wfs_dn_; }

private:
  int n_sites_;
  int n_up_;
  int n_dn_;
  arma::Mat<coeff_t> work_matrix_;
  arma::Mat<coeff_t> onebody_wfs_up_;
  arma::Mat<coeff_t> onebody_wfs_dn_;
};

} // namespace hydra
