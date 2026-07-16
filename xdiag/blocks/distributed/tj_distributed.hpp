// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_DISTRIBUTED

#include <cstdint>
#include <memory>
#include <string>

#include <xdiag/basis/basis.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/symmetries/representation_set.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

class tJDistributedIterator;

// Distributed (MPI) t-J block (no double occupancy). Number conserving in both
// spin species, no permutation symmetries. Stores only its irreps and a
// type-erased basis (a BasistJDistributed<bit_t> backend chosen from nsites).
class XDIAG_API tJDistributed {
public:
  using iterator_t = tJDistributedIterator;

  tJDistributed() = default;

  // Generic constructor
  tJDistributed(int64_t nsites, RepresentationSet const &irreps);

  // Convenience constructor
  tJDistributed(int64_t nsites, int64_t nup, int64_t ndn);

  int64_t nsites() const;
  constexpr int64_t d() const { return 3; }
  int64_t dim() const;      // global dimension
  int64_t size() const;     // locally stored dimension
  int64_t size_max() const; // largest local dimension over all ranks
  int64_t size_min() const; // smallest local dimension over all ranks
  bool isreal() const;
  int64_t index(ProductState const &pstate) const;

  bool operator==(tJDistributed const &rhs) const;
  bool operator!=(tJDistributed const &rhs) const;

  iterator_t begin() const;
  iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
};

std::ostream &operator<<(std::ostream &out, tJDistributed const &block);
std::string to_string(tJDistributed const &block);

class XDIAG_API tJDistributedIterator {
public:
  tJDistributedIterator(tJDistributed const *block, int64_t idx);
  tJDistributedIterator &operator++();
  ProductState operator*() const;
  bool operator==(tJDistributedIterator const &rhs) const;
  bool operator!=(tJDistributedIterator const &rhs) const;

private:
  std::unique_ptr<basis::BasisIterator> it_;
  int64_t idx_;
};

} // namespace xdiag
#endif
