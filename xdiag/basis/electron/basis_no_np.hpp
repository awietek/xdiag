#pragma once

#include <xdiag/common.hpp>

#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/combinatorics/subsets_index.hpp>

namespace xdiag::basis::electron {

using namespace combinatorics;

template <typename bit_tt> class BasisNoNpIterator;

template <typename bit_tt> class BasisNoNp {
public:
  using bit_t = bit_tt;
  using iterator_t = BasisNoNpIterator<bit_t>;

  BasisNoNp() = default;
  BasisNoNp(int n_sites);

  int64_t dim() const;
  int64_t size() const;
  iterator_t begin() const;
  iterator_t end() const;

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

template <typename bit_tt> class BasisNoNpIterator {
public:
  using bit_t = bit_tt;
  BasisNoNpIterator() = default;
  BasisNoNpIterator(int64_t n_sites, bool begin);
  BasisNoNpIterator &operator++();
  std::pair<bit_t, bit_t> operator*() const;
  bool operator!=(BasisNoNpIterator<bit_t> const &rhs) const;

private:
  bit_t max_;
  bit_t ups_;
  bit_t dns_;
};

} // namespace xdiag::basis::electron
