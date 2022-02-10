#pragma once

#include <hydra/common.h>

#include <hydra/combinatorics/combinations.h>

#include <hydra/indexing/lintable.h>

namespace hydra::indexing {
template <typename bit_t> class SpinhalfIndexingSz {
public:
  SpinhalfIndexingSz() = default;
  SpinhalfIndexingSz(int n_sites, int nup);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline idx_t size() const { return size_; }
  inline idx_t index(bit_t spins) const { return lintable_.index(spins); }

  inline combinatorics::Combinations<bit_t> states() const {
    return combinatorics::Combinations<bit_t>(n_sites_, n_up_);
  }

private:
  int n_sites_;
  int n_up_;
  indexing::LinTable<bit_t> lintable_;
  idx_t size_;
};

} // namespace hydra::indexing
