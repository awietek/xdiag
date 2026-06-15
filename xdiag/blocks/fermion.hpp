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

class Fermion {
public:
  using iterator_t = FermionIterator;

  XDIAG_API Fermion() = default;

  // Generic constructor
  Fermion(int64_t nsites, RepresentationSet const &irreps);

  // Convenience constructors
  XDIAG_API Fermion(int64_t nsites);
  XDIAG_API Fermion(int64_t nsites, int64_t number);
  XDIAG_API Fermion(int64_t nsites, Representation const &irrep);
  XDIAG_API Fermion(int64_t nsites, int64_t number,
                    Representation const &irrep);

  XDIAG_API int64_t nsites() const;
  XDIAG_API constexpr int64_t d() const { return 2; }
  XDIAG_API int64_t dim() const;
  XDIAG_API int64_t size() const;
  XDIAG_API bool isreal() const;

  XDIAG_API bool operator==(Fermion const &rhs) const;
  XDIAG_API bool operator!=(Fermion const &rhs) const;

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
};

XDIAG_API std::ostream &operator<<(std::ostream &out, Fermion const &block);
XDIAG_API std::string to_string(Fermion const &block);

class FermionIterator {
public:
  XDIAG_API FermionIterator(Fermion const *block, int64_t idx);
  XDIAG_API FermionIterator &operator++();
  XDIAG_API ProductState operator*() const;
  XDIAG_API bool operator==(FermionIterator const &rhs) const;
  XDIAG_API bool operator!=(FermionIterator const &rhs) const;

private:
  std::unique_ptr<basis::BasisIterator> it_;
  int64_t idx_;
};

} // namespace xdiag
