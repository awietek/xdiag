// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>

#include <xdiag/basis/basis.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/symmetries/representation_set.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

class SpinhalfIterator;

class Spinhalf {
public:
  using iterator_t = SpinhalfIterator;

  XDIAG_API Spinhalf() = default;

  // Generic constructor
  Spinhalf(int64_t sites, RepresentationSet const &irreps,
           std::string backend = "auto");

  // Convenience constructors
  XDIAG_API Spinhalf(int64_t nsites);
  XDIAG_API Spinhalf(int64_t nsites, int64_t nup);
  XDIAG_API Spinhalf(int64_t nsites, Representation const &irrep,
                     std::string backend = "auto");
  XDIAG_API Spinhalf(int64_t nsites, int64_t nup, Representation const &irrep,
                     std::string backend = "auto");

  XDIAG_API int64_t nsites() const;
  XDIAG_API int64_t dim() const;
  XDIAG_API int64_t size() const;
  XDIAG_API bool isreal() const;

  XDIAG_API bool operator==(Spinhalf const &rhs) const;
  XDIAG_API bool operator!=(Spinhalf const &rhs) const;

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  int64_t nsites_ = 0;
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
  int64_t size_ = 0;
};

XDIAG_API int64_t nsites(Spinhalf const &block);
XDIAG_API int64_t dim(Spinhalf const &block);
XDIAG_API int64_t size(Spinhalf const &block);
XDIAG_API bool isreal(Spinhalf const &block);

XDIAG_API std::ostream &operator<<(std::ostream &out, Spinhalf const &block);
XDIAG_API std::string to_string(Spinhalf const &block);

class SpinhalfIterator {
public:
  XDIAG_API SpinhalfIterator(Spinhalf const *block, int64_t idx);
  XDIAG_API SpinhalfIterator &operator++();
  XDIAG_API ProductState operator*() const;
  XDIAG_API bool operator==(SpinhalfIterator const &rhs) const;
  XDIAG_API bool operator!=(SpinhalfIterator const &rhs) const;

private:
  Spinhalf const *block_;
  int64_t idx_;
};

} // namespace xdiag
