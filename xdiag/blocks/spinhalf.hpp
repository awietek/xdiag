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

class SpinhalfIterator;

class XDIAG_API Spinhalf {
public:
  using iterator_t = SpinhalfIterator;

   Spinhalf() = default;

  // Generic constructor
  Spinhalf(int64_t sites, RepresentationSet const &irreps,
           std::string backend = "auto");

  // Convenience constructors
   Spinhalf(int64_t nsites);
   Spinhalf(int64_t nsites, int64_t nup);
   Spinhalf(int64_t nsites, Representation const &irrep,
                     std::string backend = "auto");
   Spinhalf(int64_t nsites, int64_t nup, Representation const &irrep,
                     std::string backend = "auto");

   int64_t nsites() const;
   constexpr int64_t d() const { return 2; }
   int64_t dim() const;
   int64_t size() const;
   bool isreal() const;
   int64_t index(ProductState const &pstate) const;

   bool operator==(Spinhalf const &rhs) const;
   bool operator!=(Spinhalf const &rhs) const;

   iterator_t begin() const;
   iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
};

 std::ostream &operator<<(std::ostream &out, Spinhalf const &block);
 std::string to_string(Spinhalf const &block);

class XDIAG_API SpinhalfIterator {
public:
   SpinhalfIterator(Spinhalf const *block, int64_t idx);
   SpinhalfIterator &operator++();
   ProductState operator*() const;
   bool operator==(SpinhalfIterator const &rhs) const;
   bool operator!=(SpinhalfIterator const &rhs) const;

private:
  std::unique_ptr<basis::BasisIterator> it_;
  int64_t idx_;
};

} // namespace xdiag
