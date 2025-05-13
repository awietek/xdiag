// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>

#include <xdiag/common.hpp>

#include <xdiag/basis/tj/basis_tj.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

class tJIterator;

class tJ {
public:
  using basis_t = basis::BasistJ;
  using iterator_t = tJIterator;

  XDIAG_API tJ() = default;
  XDIAG_API tJ(int64_t nsites, int64_t nup, int64_t ndn,
               std::string backend = "auto");
  XDIAG_API tJ(int64_t nsites, int64_t nup, int64_t ndn,
               Representation const &irrep, std::string backend = "auto");

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;
  XDIAG_API int64_t index(ProductState const &pstate) const;
  XDIAG_API int64_t dim() const;
  XDIAG_API int64_t size() const;

  XDIAG_API bool operator==(tJ const &rhs) const;
  XDIAG_API bool operator!=(tJ const &rhs) const;

  XDIAG_API int64_t nsites() const;
  XDIAG_API bool isreal() const;

  std::string backend() const;
  std::optional<int64_t> nup() const;
  std::optional<int64_t> ndn() const;
  std::optional<Representation> const &irrep() const;
  basis_t const &basis() const;

private:
  int64_t nsites_;
  std::string backend_;
  std::optional<int64_t> nup_;
  std::optional<int64_t> ndn_;
  std::optional<Representation> irrep_;
  std::shared_ptr<basis_t> basis_;
  int64_t size_;
};

XDIAG_API int64_t index(tJ const &block, ProductState const &pstate);
XDIAG_API int64_t nsites(tJ const &block);
XDIAG_API int64_t dim(tJ const &block);
XDIAG_API int64_t size(tJ const &block);
XDIAG_API bool isreal(tJ const &block);
XDIAG_API std::ostream &operator<<(std::ostream &out, tJ const &block);
XDIAG_API std::string to_string(tJ const &block);

class tJIterator {
public:
  tJIterator(tJ const &block, bool begin);
  XDIAG_API tJIterator &operator++();
  XDIAG_API ProductState const &operator*() const;
  XDIAG_API bool operator!=(tJIterator const &rhs) const;

private:
  int64_t nsites_;
  mutable ProductState pstate_;
  basis::BasistJIterator it_;
};

} // namespace xdiag
