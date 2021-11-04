#pragma once

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>

namespace hydra {

template <class bit_t = std_bit_t> class tJ {
public:
  tJ() = default;
  tJ(int n_sites, int nup, int ndn);

  int n_sites() const { return n_sites_; }
  int n_up() const { return n_up_; }
  int n_dn() const { return n_dn_; }

  bool charge_conserved() const { return charge_conserved_; }
  bool sz_conserved() const { return sz_conserved_; }

  idx_t size() const { return size_; }

  bool operator==(tJ const &rhs) const;
  bool operator!=(tJ const &rhs) const;

  // private:
  int n_sites_;
  bool charge_conserved_;
  int charge_;
  bool sz_conserved_;
  int sz_;
  int n_up_;
  int n_dn_;
  int n_holes_;

  idx_t size_holes_;
  idx_t size_spins_;
  idx_t size_;
  
  indexing::LinTable<bit_t> lintable_holes_;
  indexing::LinTable<bit_t> lintable_spins_;

};

} // namespace hydra
