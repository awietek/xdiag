#pragma once

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>

namespace hydra::indexing {

template <typename bit_t> class tJIndexing {
public:
  tJIndexing() = default;
  tJIndexing(int n_sites, int n_up, int n_dn);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }
  inline int n_spins() const { return n_spins_; }
  inline int n_holes() const { return n_holes_; }

  inline idx_t size_spins() const { return size_spins_; }
  inline idx_t size_holes() const { return size_holes_; }

  inline idx_t size() const { return size_; }
  inline idx_t index_holes(bit_t holes) const {
    return lintable_holes_.index(holes);
  }
  inline idx_t index_spins(bit_t spins) const {
    return lintable_spins_.index(spins);
  }

  // For testing purposes
  inline idx_t index(bit_t ups, bit_t dns) const {
    bit_t sitesmask = ((bit_t)1 << n_sites_) - 1;
    bit_t not_holes = ups | dns;
    bit_t holes = (~not_holes) & sitesmask;
    bit_t spins = bitops::extract(ups, not_holes);
    return index_holes(holes) * size_spins_ + index_spins(spins);
  }



private:
  int n_sites_;
  int n_up_;
  int n_dn_;
  int n_spins_;
  int n_holes_;

  idx_t size_holes_;
  idx_t size_spins_;
  idx_t size_;

  indexing::LinTable<bit_t> lintable_holes_;
  indexing::LinTable<bit_t> lintable_spins_;
};

} // namespace hydra::indexing
