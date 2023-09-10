#pragma once

#include <hydra/common.h>

#include <hydra/combinatorics/subsets.h>
#include <hydra/combinatorics/subsets_index.h>

namespace hydra::basis::electron {

using namespace combinatorics;

template <typename bit_t> class BasisNoNp {
public:
  using bit_type = bit_t;
  
  BasisNoNp() = default;
  BasisNoNp(int n_sites);

  int64_t size() const;
  inline int64_t index(bit_t ups, bit_t dns) const {
    return index_ups(ups) * size_dns_ + index_dns(dns);
  }

  int n_sites() const;
  static constexpr bool np_conserved() { return false; }

private:
  int n_sites_;
  int64_t size_ups_;
  int64_t size_dns_;
  int64_t size_;

public:
  int64_t size_ups() const;
  int64_t size_dns() const;
  inline int64_t index_ups(bit_t ups) const { return (int64_t)ups; }
  inline int64_t index_dns(bit_t dns) const { return (int64_t)dns; }

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
};

} // namespace hydra::basis::electron
