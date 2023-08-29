#pragma once

#include <hydra/common.h>

#include <hydra/combinatorics/subsets_index.h>

namespace hydra::basis::spinhalf {

template <typename bit_t> class BasisNoSz {
public:
  using iterator_t = combinatorics::SubsetsIterator<bit_t>;
  using bit_type = bit_t;
  
  BasisNoSz() = default;
  BasisNoSz(int64_t n_sites);

  iterator_t begin() const;
  iterator_t end() const;
  int64_t size() const;
  inline int64_t index(bit_t spins) const { return (int64_t)spins; }
  inline bit_t state(int64_t index) const { return (bit_t)index; }

  int64_t n_sites() const;
  static constexpr bool sz_conserved() { return false; }

  bool operator==(BasisNoSz const &rhs) const;
  bool operator!=(BasisNoSz const &rhs) const;

private:
  int64_t n_sites_;
  int64_t size_;
  iterator_t begin_, end_;
};

} // namespace hydra::basis::spinhalf
