#pragma once

#include <lila/all.h>

#include <hydra/common.h>

namespace hydra {

template <class bit_t, class coeff_t> class GPWFSpinhalf {
public:
  GPWFSpinhalf(int n_sites, lila::Matrix<coeff_t> const &onebody_wfs, int n_up = -1);
  int n_sites() const { return n_sites_; }
  int n_up() const { return n_up_; }
  int n_dn() const { return n_dn_; }
  coeff_t coefficient(bit_t const &state) const;

private:
  int n_sites_;
  int n_up_;
  int n_dn_;
  lila::Matrix<coeff_t> work_matrix_;
  lila::Matrix<coeff_t> onebody_wfs_up_;
  lila::Matrix<coeff_t> onebody_wfs_dn_;
};

} // namespace hydra
