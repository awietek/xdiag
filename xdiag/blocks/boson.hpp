// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <memory>

#include <xdiag/basis/basis.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/symmetries/representation_set.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

class BosonIterator;

class XDIAG_API Boson {
public:
  using iterator_t = BosonIterator;

  Boson() = default;

  // Generic constructor
  Boson(int64_t sites, int64_t d, RepresentationSet const &irreps);

  // Convenience constructors
  Boson(int64_t sites, int64_t d);
  Boson(int64_t sites, int64_t d, int64_t number);
  Boson(int64_t sites, int64_t d, Representation const &irrep);
  Boson(int64_t sites, int64_t d, int64_t number, Representation const &irrep);

  int64_t nsites() const;
  int64_t d() const;
  int64_t dim() const;
  int64_t size() const;
  bool isreal() const;
  int64_t index(ProductState const &pstate) const;

  bool operator==(Boson const &rhs) const;
  bool operator!=(Boson const &rhs) const;

  iterator_t begin() const;
  iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  int64_t d_;
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
};

XDIAG_API std::ostream &operator<<(std::ostream &out, Boson const &block);
XDIAG_API std::string to_string(Boson const &block);

class XDIAG_API BosonIterator {
public:
  BosonIterator(Boson const *block, int64_t idx);
  BosonIterator &operator++();
  ProductState operator*() const;
  bool operator==(BosonIterator const &rhs) const;
  bool operator!=(BosonIterator const &rhs) const;

private:
  std::unique_ptr<basis::BasisIterator> it_;
  int64_t idx_;
};

using Spin = Boson;

} // namespace xdiag
