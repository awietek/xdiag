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

class SpinhalfDistributedIterator;

// Distributed (MPI) spin-1/2 block. Number (Sz) conserving only, no permutation
// symmetries. Like the serial blocks it stores only its irreps and a
// type-erased basis; the basis is one of the BasisSpinhalfDistributed<bit_t>
// backends, with the bit width chosen automatically from nsites.
class XDIAG_API SpinhalfDistributed {
public:
  using iterator_t = SpinhalfDistributedIterator;

  SpinhalfDistributed() = default;

  // Generic constructor
  SpinhalfDistributed(int64_t nsites, RepresentationSet const &irreps);

  // Convenience constructor
  SpinhalfDistributed(int64_t nsites, int64_t nup);

  int64_t nsites() const;
  constexpr int64_t d() const { return 2; }
  int64_t dim() const;      // global dimension
  int64_t size() const;     // locally stored dimension
  int64_t size_max() const; // largest local dimension over all ranks
  int64_t size_min() const; // smallest local dimension over all ranks
  bool isreal() const;
  int64_t index(ProductState const &pstate) const;

  bool operator==(SpinhalfDistributed const &rhs) const;
  bool operator!=(SpinhalfDistributed const &rhs) const;

  iterator_t begin() const;
  iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
};

std::ostream &operator<<(std::ostream &out, SpinhalfDistributed const &block);
std::string to_string(SpinhalfDistributed const &block);

class XDIAG_API SpinhalfDistributedIterator {
public:
  SpinhalfDistributedIterator(SpinhalfDistributed const *block, int64_t idx);
  SpinhalfDistributedIterator &operator++();
  ProductState operator*() const;
  bool operator==(SpinhalfDistributedIterator const &rhs) const;
  bool operator!=(SpinhalfDistributedIterator const &rhs) const;

private:
  std::unique_ptr<basis::BasisIterator> it_;
  int64_t idx_;
};

} // namespace xdiag
#endif
