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

class Electron {
public:
  using iterator_t = ElectronIterator;

  XDIAG_API Electron() = default;

  // Generic constructor
  Electron(int64_t sites, RepresentationSet const &irreps);

  // Convenience constructors
  XDIAG_API Electron(int64_t nsites);
  XDIAG_API Electron(int64_t nsites, int64_t nup, int64_t ndn);
  XDIAG_API Electron(int64_t nsites, Representation const &irrep);
  XDIAG_API Electron(int64_t nsites, int64_t nup, int64_t ndn,
                     Representation const &irrep);

  XDIAG_API int64_t nsites() const;
  XDIAG_API constexpr int64_t d() const { return 4; }
  XDIAG_API int64_t dim() const;
  XDIAG_API int64_t size() const;
  XDIAG_API bool isreal() const;

  XDIAG_API bool operator==(Electron const &rhs) const;
  XDIAG_API bool operator!=(Electron const &rhs) const;

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;

  RepresentationSet irreps() const;
  std::shared_ptr<basis::Basis> const &basis() const;

private:
  RepresentationSet irreps_;
  std::shared_ptr<basis::Basis> basis_;
};

XDIAG_API std::ostream &operator<<(std::ostream &out, Electron const &block);
XDIAG_API std::string to_string(Electron const &block);

class ElectronIterator {
public:
  XDIAG_API ElectronIterator(Electron const *block, int64_t idx);
  XDIAG_API ElectronIterator &operator++();
  XDIAG_API ProductState operator*() const;
  XDIAG_API bool operator==(ElectronIterator const &rhs) const;
  XDIAG_API bool operator!=(ElectronIterator const &rhs) const;

private:
  std::unique_ptr<basis::BasisIterator> it_;
  int64_t idx_;
};

} // namespace xdiag
