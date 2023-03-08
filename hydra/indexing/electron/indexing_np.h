#pragma once

#include <hydra/common.h>

#include <hydra/combinatorics/combinations.h>

#include <hydra/indexing/lin_table.h>

namespace hydra::indexing::electron {

template <typename bit_t> class IndexingNp {
public:
  IndexingNp() = default;
  IndexingNp(int n_sites, int n_up, int n_dn);

  int n_sites() const;
  int n_up() const;
  int n_dn() const;
  idx_t size_ups() const;
  idx_t size_dns() const;
  idx_t size() const;

  inline idx_t index_ups(bit_t ups) const { return lintable_ups_.index(ups); }
  inline idx_t index_dns(bit_t dns) const { return lintable_dns_.index(dns); }
  inline idx_t index(bit_t ups, bit_t dns) const {
    return index_ups(ups) * size_dns_ + index_dns(dns);
  }

  combinatorics::Combinations<bit_t> states_ups() const;
  combinatorics::Combinations<bit_t> states_dns() const;

private:
  int n_sites_;
  int n_up_;
  int n_dn_;

  idx_t size_ups_;
  idx_t size_dns_;
  idx_t size_;

  indexing::LinTable<bit_t> lintable_ups_;
  indexing::LinTable<bit_t> lintable_dns_;
};

} // namespace hydra::indexing::electron
