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

class ElectronDistributedIterator;

// Distributed (MPI) electron block (double occupancy allowed). Number
// conserving in both spin species, no permutation symmetries. Stores only its
// irreps and a type-erased basis (a BasisElectronDistributed<bit_t> backend
// from nsites).
class XDIAG_API ElectronDistributed {
public:
  using iterator_t = ElectronDistributedIterator;

  ElectronDistributed() = default;

  // Generic constructor
  ElectronDistributed(int64_t nsites, RepresentationSet const &irreps);

  // Convenience constructor
  ElectronDistributed(int64_t nsites, int64_t nup, int64_t ndn);

  int64_t nsites() const;
  constexpr int64_t d() const { return 4; }
  int64_t dim() const;      // global dimension
  int64_t size() const;     // locally stored dimension
  int64_t size_max() const; // largest local dimension over all ranks
  int64_t size_min() const; // smallest local dimension over all ranks
  bool isreal() const;
  int64_t index(ProductState const &pstate) const;

  bool operator==(ElectronDistributed const &rhs) const;
  bool operator!=(ElectronDistributed const &rhs) const;

  iterator_t begin() const;
  iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
};

std::ostream &operator<<(std::ostream &out, ElectronDistributed const &block);
std::string to_string(ElectronDistributed const &block);

class XDIAG_API ElectronDistributedIterator {
public:
  ElectronDistributedIterator(ElectronDistributed const *block, int64_t idx);
  ElectronDistributedIterator &operator++();
  ProductState operator*() const;
  bool operator==(ElectronDistributedIterator const &rhs) const;
  bool operator!=(ElectronDistributedIterator const &rhs) const;

private:
  std::unique_ptr<basis::BasisIterator> it_;
  int64_t idx_;
};

} // namespace xdiag
#endif
