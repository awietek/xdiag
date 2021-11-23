#pragma once

#include <hydra/common.h>
#include <hydra/indexing/tj/tj_indexing.h>

namespace hydra {

template <class bit_t> class tJ {
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

  indexing::tJIndexing<bit_t> const &indexing() const { return indexing_; }

private:
  int n_sites_;
  bool charge_conserved_;
  int charge_;
  bool sz_conserved_;
  int sz_;
  int n_up_;
  int n_dn_;

  indexing::tJIndexing<bit_t> indexing_;
  idx_t size_;
};

} // namespace hydra
