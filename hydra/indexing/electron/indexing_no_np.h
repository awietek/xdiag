#pragma once

#include <hydra/common.h>

#include <hydra/combinatorics/subsets.h>
#include <hydra/combinatorics/subsets_index.h>

namespace hydra::indexing::electron {

using namespace combinatorics;

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

  Subsets<bit_t> states_ups() const;
  Subsets<bit_t> states_dns() const;
  SubsetsIndex<bit_t> states_indices_ups() const;
  SubsetsIndex<bit_t> states_indices_dns() const;

#ifdef _OPENMP
  SubsetsThread<bit_t> states_ups_thread() const;
  SubsetsThread<bit_t> states_dns_thread() const;
  SubsetsIndexThread<bit_t> states_indices_ups_thread() const;
  SubsetsIndexThread<bit_t> states_indices_dns_thread() const;
#endif

private:
  int n_sites_;
  idx_t size_ups_;
  idx_t size_dns_;
  idx_t size_;
};

} // namespace hydra::indexing::electron
