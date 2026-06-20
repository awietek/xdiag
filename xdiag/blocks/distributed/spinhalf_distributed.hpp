// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

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
// symmetries. Like the serial blocks it stores only its irreps and a type-erased
// basis; the basis is one of the BasisSpinhalfDistributed<bit_t> backends, with
// the bit width chosen automatically from nsites.
class SpinhalfDistributed {
public:
  using iterator_t = SpinhalfDistributedIterator;

  XDIAG_API SpinhalfDistributed() = default;

  // Generic constructor
  SpinhalfDistributed(int64_t nsites, RepresentationSet const &irreps);

  // Convenience constructor
  XDIAG_API SpinhalfDistributed(int64_t nsites, int64_t nup);

  XDIAG_API int64_t nsites() const;
  XDIAG_API constexpr int64_t d() const { return 2; }
  XDIAG_API int64_t dim() const;       // global dimension
  XDIAG_API int64_t size() const;      // locally stored dimension
  XDIAG_API int64_t size_max() const;  // largest local dimension over all ranks
  XDIAG_API int64_t size_min() const;  // smallest local dimension over all ranks
  XDIAG_API bool isreal() const;

  XDIAG_API bool operator==(SpinhalfDistributed const &rhs) const;
  XDIAG_API bool operator!=(SpinhalfDistributed const &rhs) const;

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
};

XDIAG_API std::ostream &operator<<(std::ostream &out,
                                   SpinhalfDistributed const &block);
XDIAG_API std::string to_string(SpinhalfDistributed const &block);

class SpinhalfDistributedIterator {
public:
  XDIAG_API SpinhalfDistributedIterator(SpinhalfDistributed const *block,
                                        int64_t idx);
  XDIAG_API SpinhalfDistributedIterator &operator++();
  XDIAG_API ProductState operator*() const;
  XDIAG_API bool operator==(SpinhalfDistributedIterator const &rhs) const;
  XDIAG_API bool operator!=(SpinhalfDistributedIterator const &rhs) const;

private:
  std::unique_ptr<basis::BasisIterator> it_;
  int64_t idx_;
};

} // namespace xdiag
