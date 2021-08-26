#pragma once

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>

namespace hydra {

template <class bit_t = std_bit_t> class tJ {
public:
  tJ() = default;
  tJ(int n_sites, int nup, int ndn);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }

  inline bool charge_conserved() const { return charge_conserved_; }
  inline bool sz_conserved() const { return sz_conserved_; }

  inline idx_t size() const { return size_; }
  idx_t index(bit_t up, bit_t dn) const;

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

  idx_t size_up_;
  idx_t size_holes_;
  idx_t size_;

  std::vector<std::pair<idx_t, idx_t>> dn_limits_for_up_;
  std::vector<bit_t> dns_;
  LinTable<bit_t> lintable_up_;
};

} // namespace hydra
