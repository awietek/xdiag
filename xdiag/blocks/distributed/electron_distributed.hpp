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

class ElectronDistributedIterator;

// Distributed (MPI) electron block (double occupancy allowed). Number conserving
// in both spin species, no permutation symmetries. Stores only its irreps and a
// type-erased basis (a BasisElectronDistributed<bit_t> backend from nsites).
class ElectronDistributed {
public:
  using iterator_t = ElectronDistributedIterator;

  XDIAG_API ElectronDistributed() = default;

  // Generic constructor
  ElectronDistributed(int64_t nsites, RepresentationSet const &irreps);

  // Convenience constructor
  XDIAG_API ElectronDistributed(int64_t nsites, int64_t nup, int64_t ndn);

  XDIAG_API int64_t nsites() const;
  XDIAG_API constexpr int64_t d() const { return 4; }
  XDIAG_API int64_t dim() const;       // global dimension
  XDIAG_API int64_t size() const;      // locally stored dimension
  XDIAG_API int64_t size_max() const;  // largest local dimension over all ranks
  XDIAG_API int64_t size_min() const;  // smallest local dimension over all ranks
  XDIAG_API bool isreal() const;
  XDIAG_API int64_t index(ProductState const &pstate) const;

  XDIAG_API bool operator==(ElectronDistributed const &rhs) const;
  XDIAG_API bool operator!=(ElectronDistributed const &rhs) const;

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
};

XDIAG_API std::ostream &operator<<(std::ostream &out,
                                   ElectronDistributed const &block);
XDIAG_API std::string to_string(ElectronDistributed const &block);

class ElectronDistributedIterator {
public:
  XDIAG_API ElectronDistributedIterator(ElectronDistributed const *block,
                                        int64_t idx);
  XDIAG_API ElectronDistributedIterator &operator++();
  XDIAG_API ProductState operator*() const;
  XDIAG_API bool operator==(ElectronDistributedIterator const &rhs) const;
  XDIAG_API bool operator!=(ElectronDistributedIterator const &rhs) const;

private:
  std::unique_ptr<basis::BasisIterator> it_;
  int64_t idx_;
};

} // namespace xdiag
