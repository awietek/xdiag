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

class Boson {
public:
  using iterator_t = BosonIterator;

  XDIAG_API Boson() = default;

  // Generic constructor
  Boson(int64_t sites, int64_t d, RepresentationSet const &irreps);

  // Convenience constructors
  XDIAG_API Boson(int64_t sites, int64_t d);
  XDIAG_API Boson(int64_t sites, int64_t d, int64_t number);
  XDIAG_API Boson(int64_t sites, int64_t d, Representation const &irrep);
  XDIAG_API Boson(int64_t sites, int64_t d, int64_t number,
                  Representation const &irrep);

  XDIAG_API int64_t nsites() const;
  XDIAG_API int64_t d() const;
  XDIAG_API int64_t dim() const;
  XDIAG_API int64_t size() const;
  XDIAG_API bool isreal() const;
  XDIAG_API int64_t index(ProductState const &pstate) const;

  XDIAG_API bool operator==(Boson const &rhs) const;
  XDIAG_API bool operator!=(Boson const &rhs) const;

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  int64_t d_;
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
};

XDIAG_API std::ostream &operator<<(std::ostream &out, Boson const &block);
XDIAG_API std::string to_string(Boson const &block);

class BosonIterator {
public:
  XDIAG_API BosonIterator(Boson const *block, int64_t idx);
  XDIAG_API BosonIterator &operator++();
  XDIAG_API ProductState operator*() const;
  XDIAG_API bool operator==(BosonIterator const &rhs) const;
  XDIAG_API bool operator!=(BosonIterator const &rhs) const;

private:
  std::unique_ptr<basis::BasisIterator> it_;
  int64_t idx_;
};

using Spin = Boson;

} // namespace xdiag
