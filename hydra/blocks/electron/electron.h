#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/blocks/blocks.h>
#include <hydra/indexing/electron/electron_indexing.h>

namespace hydra {

template <class bit_t> class Electron {
public:
  Electron() = default;
  Electron(int n_sites, int nup, int ndn);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }

  inline bool charge_conserved() const { return charge_conserved_; }
  inline bool sz_conserved() const { return sz_conserved_; }
  inline idx_t size() const { return size_; }

  bool operator==(Electron const &rhs) const;
  bool operator!=(Electron const &rhs) const;
  indexing::ElectronIndexing<bit_t> const &indexing() const { return indexing_; }

private:
  int n_sites_;

  bool charge_conserved_;
  int charge_;
  bool sz_conserved_;
  int sz_;
  int n_up_;
  int n_dn_;

  indexing::ElectronIndexing<bit_t> indexing_;

  idx_t size_;
};

} // namespace hydra
