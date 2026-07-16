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

class tJIterator;

// Spinful t-J block (local dimension d = 3: empty, up, dn). Like Electron but
// with the no-double-occupancy constraint; the dn sector is stored compressed
// into the non-up sites (see basis/basis_tj.hpp). Mirrors the Electron block's
// RepresentationSet-based interface.
class XDIAG_API tJ {
public:
  using iterator_t = tJIterator;

  tJ() = default;

  // Generic constructor
  tJ(int64_t nsites, RepresentationSet const &irreps);

  // Convenience constructors
  tJ(int64_t nsites);
  tJ(int64_t nsites, int64_t nup, int64_t ndn);
  tJ(int64_t nsites, Representation const &irrep);
  tJ(int64_t nsites, int64_t nup, int64_t ndn, Representation const &irrep);

  int64_t nsites() const;
  constexpr int64_t d() const { return 3; }
  int64_t dim() const;
  int64_t size() const;
  bool isreal() const;
  int64_t index(ProductState const &pstate) const;

  bool operator==(tJ const &rhs) const;
  bool operator!=(tJ const &rhs) const;

  iterator_t begin() const;
  iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
};

XDIAG_API std::ostream &operator<<(std::ostream &out, tJ const &block);
XDIAG_API std::string to_string(tJ const &block);

class XDIAG_API tJIterator {
public:
  tJIterator(tJ const *block, int64_t idx);
  tJIterator &operator++();
  ProductState operator*() const;
  bool operator==(tJIterator const &rhs) const;
  bool operator!=(tJIterator const &rhs) const;

private:
  std::unique_ptr<basis::BasisIterator> it_;
  int64_t idx_;
};

} // namespace xdiag
