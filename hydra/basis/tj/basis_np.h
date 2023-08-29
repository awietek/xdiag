#pragma once

#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/lin_table.h>
#include <hydra/common.h>

namespace hydra::basis::tj {

template <typename bit_t> class BasisNp {
public:
  using bit_type = bit_t;

  BasisNp() = default;
  BasisNp(int n_sites, int n_up, int n_dn);

  idx_t size() const;
  inline idx_t index(bit_t ups, bit_t dns) const {
    bit_t dncs = bits::extract<bit_t>(dns, (~ups) & sitesmask_);
    return index_ups(ups) * size_dncs_ + index_dncs(dncs);
  }

  int n_sites() const;
  int n_up() const;
  int n_dn() const;
  static constexpr bool np_conserved() { return true; }

private:
  int n_sites_;
  int n_up_;
  int n_dn_;

  idx_t size_ups_;
  idx_t size_dncs_;
  idx_t size_;
  bit_t sitesmask_;

  combinatorics::LinTable<bit_t> lintable_ups_;
  combinatorics::LinTable<bit_t> lintable_dncs_;

public:
  idx_t size_ups() const;
  idx_t size_dncs() const;
  idx_t ups_offset(idx_t idx_ups) const;

  combinatorics::Combinations<bit_t> states_ups() const;
  combinatorics::Combinations<bit_t> states_dncs(bit_t ups) const;
  combinatorics::CombinationsIndex<bit_t> states_indices_ups() const;
  combinatorics::CombinationsIndex<bit_t> states_indices_dncs(bit_t ups) const;

#ifdef _OPENMP
  combinatorics::CombinationsThread<bit_t> states_ups_thread() const;
  combinatorics::CombinationsThread<bit_t> states_dncs_thread(bit_t ups) const;
  combinatorics::CombinationsIndexThread<bit_t>
  states_indices_ups_thread() const;
  combinatorics::CombinationsIndexThread<bit_t>
  states_indices_dncs_thread(bit_t ups) const;
#endif

  inline idx_t index_ups(bit_t ups) const { return lintable_ups_.index(ups); }
  inline idx_t index_dncs(bit_t dncs) const {
    return lintable_dncs_.index(dncs);
  }
};

} // namespace hydra::basis::tj
