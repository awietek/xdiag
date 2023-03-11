#pragma once

#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/indexing/lin_table.h>

namespace hydra::indexing::tj {

using namespace combinatorics;

template <typename bit_t> class IndexingNp {
public:
  IndexingNp() = default;
  IndexingNp(int n_sites, int n_up, int n_dn);

  int n_sites() const;
  int n_up() const;
  int n_dn() const;
  idx_t size_ups() const;
  idx_t size_dncs() const;
  idx_t size() const;

  inline idx_t index_ups(bit_t ups) const { return lintable_ups_.index(ups); }
  inline idx_t index_dncs(bit_t dncs) const {
    return lintable_dncs_.index(dncs);
  }
  inline idx_t index(bit_t ups, bit_t dns) const {
    bit_t dncs = bitops::extract<bit_t>(dns, (~ups) & sitesmask_);
    return index_ups(ups) * size_dncs_ + index_dncs(dncs);
  }
  

  Combinations<bit_t> states_ups() const;
  Combinations<bit_t> states_dncs(bit_t ups) const;
  CombinationsIndex<bit_t> states_indices_ups() const;
  CombinationsIndex<bit_t> states_indices_dncs(bit_t ups) const;

#ifdef _OPENMP
  CombinationsThread<bit_t> states_ups_thread() const;
  CombinationsThread<bit_t> states_dncs_thread(bit_t ups) const;
  CombinationsIndexThread<bit_t> states_indices_ups_thread() const;
  CombinationsIndexThread<bit_t> states_indices_dncs_thread(bit_t ups) const;
#endif

private:
  int n_sites_;
  int n_up_;
  int n_dn_;

  idx_t size_ups_;
  idx_t size_dncs_;
  idx_t size_;

  bit_t sitesmask_;

  indexing::LinTable<bit_t> lintable_ups_;
  indexing::LinTable<bit_t> lintable_dncs_;
};

} // namespace hydra::indexing::tj
