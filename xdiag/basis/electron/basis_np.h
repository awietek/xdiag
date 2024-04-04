#pragma once

#include <xdiag/common.h>

#include <xdiag/combinatorics/combinations.h>
#include <xdiag/combinatorics/combinations_index.h>
#include <xdiag/combinatorics/lin_table.h>

namespace xdiag::basis::electron {

using namespace combinatorics;

template <typename bit_t> class BasisNp {
public:
  using bit_type = bit_t;
  
  BasisNp() = default;
  BasisNp(int n_sites, int n_up, int n_dn);

  int64_t size() const;
  int64_t dim() const;

  inline int64_t index(bit_t ups, bit_t dns) const {
    return index_ups(ups) * size_dns_ + index_dns(dns);
  }

  int n_sites() const;
  int n_up() const;
  int n_dn() const;
  static constexpr bool np_conserved() { return true; }

private:
  int n_sites_;
  int n_up_;
  int n_dn_;

  int64_t size_ups_;
  int64_t size_dns_;
  int64_t size_;

  combinatorics::LinTable<bit_t> lintable_ups_;
  combinatorics::LinTable<bit_t> lintable_dns_;

public:
  Combinations<bit_t> states_ups() const;
  Combinations<bit_t> states_dns() const;
  CombinationsIndex<bit_t> states_indices_ups() const;
  CombinationsIndex<bit_t> states_indices_dns() const;

#ifdef _OPENMP
  CombinationsThread<bit_t> states_ups_thread() const;
  CombinationsThread<bit_t> states_dns_thread() const;
  CombinationsIndexThread<bit_t> states_indices_ups_thread() const;
  CombinationsIndexThread<bit_t> states_indices_dns_thread() const;
#endif

  int64_t size_ups() const;
  int64_t size_dns() const;
  inline int64_t index_ups(bit_t ups) const { return lintable_ups_.index(ups); }
  inline int64_t index_dns(bit_t dns) const { return lintable_dns_.index(dns); }
};

} // namespace xdiag::basis::electron
