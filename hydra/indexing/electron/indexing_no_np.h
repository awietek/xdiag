#pragma once

#include <hydra/common.h>

#include <hydra/combinatorics/subsets.h>

namespace hydra::indexing::electron {

template <typename bit_t> class IndexingNoNp {
public:
  IndexingNoNp() = default;
  IndexingNoNp(int n_sites);

  int n_sites() const;
  idx_t size_ups() const;
  idx_t size_dns() const;
  idx_t size() const;

  inline idx_t index_ups(bit_t ups) const { return (idx_t)ups; }
  inline idx_t index_dns(bit_t dns) const { return (idx_t)dns; }
  inline idx_t index(bit_t ups, bit_t dns) const {
    return index_ups(ups) * size_dns_ + index_dns(dns);
  }

  combinatorics::Subsets<bit_t> states_ups() const;
  combinatorics::Subsets<bit_t> states_dns() const;

private:
  int n_sites_;
  idx_t size_ups_;
  idx_t size_dns_;
  idx_t size_;
};

} // namespace hydra::indexing::electron
