// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
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

class ElectronIterator;

class XDIAG_API Electron {
public:
  using iterator_t = ElectronIterator;

  Electron() = default;

  // Generic constructor
  Electron(int64_t sites, RepresentationSet const &irreps);

  // Convenience constructors
  Electron(int64_t nsites);
  Electron(int64_t nsites, int64_t nup, int64_t ndn);
  Electron(int64_t nsites, Representation const &irrep);
  Electron(int64_t nsites, int64_t nup, int64_t ndn,
           Representation const &irrep);

  int64_t nsites() const;
  constexpr int64_t d() const { return 4; }
  int64_t dim() const;
  int64_t size() const;
  bool isreal() const;
  int64_t index(ProductState const &pstate) const;

  bool operator==(Electron const &rhs) const;
  bool operator!=(Electron const &rhs) const;

  iterator_t begin() const;
  iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
};

XDIAG_API std::ostream &operator<<(std::ostream &out, Electron const &block);
XDIAG_API std::string to_string(Electron const &block);

class XDIAG_API ElectronIterator {
public:
  ElectronIterator(Electron const *block, int64_t idx);
  ElectronIterator &operator++();
  ProductState operator*() const;
  bool operator==(ElectronIterator const &rhs) const;
  bool operator!=(ElectronIterator const &rhs) const;

private:
  std::unique_ptr<basis::BasisIterator> it_;
  int64_t idx_;
};

} // namespace xdiag
