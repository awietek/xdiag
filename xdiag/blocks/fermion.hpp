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

class FermionIterator;

class XDIAG_API Fermion {
public:
  using iterator_t = FermionIterator;

  Fermion() = default;

  // Generic constructor
  Fermion(int64_t nsites, RepresentationSet const &irreps);

  // Convenience constructors
  Fermion(int64_t nsites);
  Fermion(int64_t nsites, int64_t number);
  Fermion(int64_t nsites, Representation const &irrep);
  Fermion(int64_t nsites, int64_t number, Representation const &irrep);

  int64_t nsites() const;
  constexpr int64_t d() const { return 2; }
  int64_t dim() const;
  int64_t size() const;
  bool isreal() const;
  int64_t index(ProductState const &pstate) const;

  bool operator==(Fermion const &rhs) const;
  bool operator!=(Fermion const &rhs) const;

  iterator_t begin() const;
  iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
};

std::ostream &operator<<(std::ostream &out, Fermion const &block);
std::string to_string(Fermion const &block);

class XDIAG_API FermionIterator {
public:
  FermionIterator(Fermion const *block, int64_t idx);
  FermionIterator &operator++();
  ProductState operator*() const;
  bool operator==(FermionIterator const &rhs) const;
  bool operator!=(FermionIterator const &rhs) const;

private:
  std::unique_ptr<basis::BasisIterator> it_;
  int64_t idx_;
};

} // namespace xdiag
