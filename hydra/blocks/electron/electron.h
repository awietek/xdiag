#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/blocks/blocks.h>
#include <hydra/indexing/lintable.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t> class Electron {
public:
  Electron() = default;
  Electron(int n_sites, int nup, int ndn);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }

  inline indexing::LinTable<bit_t> const &lintable_up() const {
    return lintable_up_;
  }
  inline indexing::LinTable<bit_t> const &lintable_dn() const {
    return lintable_dn_;
  }

  inline bool charge_conserved() const { return charge_conserved_; }
  inline bool sz_conserved() const { return sz_conserved_; }

  inline idx_t size_up() const { return size_up_; }
  inline idx_t size_dn() const { return size_dn_; }
  inline idx_t size() const { return size_; }

  inline idx_t index_up(bit_t up) const { return lintable_up_.index(up); }
  inline idx_t index_dn(bit_t dn) const { return lintable_dn_.index(dn); }
  inline idx_t index(bit_t up, bit_t dn) const {
    return index_up(up) * size_dn_ + index_dn(dn);
  }

  bool operator==(Electron const &rhs) const;
  bool operator!=(Electron const &rhs) const;

private:
  int n_sites_;

  bool charge_conserved_;
  int charge_;
  bool sz_conserved_;
  int sz_;
  int n_up_;
  int n_dn_;

  indexing::LinTable<bit_t> lintable_up_;
  indexing::LinTable<bit_t> lintable_dn_;

  idx_t size_up_;
  idx_t size_dn_;
  idx_t size_;
};

} // namespace hydra
