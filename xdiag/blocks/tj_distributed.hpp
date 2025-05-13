// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_USE_MPI
#include <optional>

#include <xdiag/basis/tj_distributed/basis_tj_distributed.hpp>
#include <xdiag/common.hpp>
#include <xdiag/states/product_state.hpp>

namespace xdiag {

class tJDistributedIterator;

class tJDistributed {
public:
  using basis_t = basis::BasistJDistributed;
  using iterator_t = tJDistributedIterator;

  XDIAG_API tJDistributed() = default;
  XDIAG_API tJDistributed(int64_t nsites, int64_t nup, int64_t ndn,
                          std::string backend = "auto");

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;
  XDIAG_API int64_t index(ProductState const &pstate) const;
  XDIAG_API int64_t dim() const;
  XDIAG_API int64_t size() const;
  XDIAG_API int64_t size_max() const;
  XDIAG_API int64_t size_min() const;

  XDIAG_API bool operator==(tJDistributed const &rhs) const;
  XDIAG_API bool operator!=(tJDistributed const &rhs) const;

  XDIAG_API int64_t nsites() const;
  XDIAG_API bool isreal() const;

  std::string backend() const;
  std::optional<int64_t> nup() const;
  std::optional<int64_t> ndn() const;
  basis_t const &basis() const;

private:
  int64_t nsites_;
  std::string backend_;
  std::optional<int64_t> nup_;
  std::optional<int64_t> ndn_;
  std::shared_ptr<basis_t> basis_;
  int64_t dim_;
  int64_t size_;
};

XDIAG_API int64_t index(tJDistributed const &block, ProductState const &pstate);
XDIAG_API int64_t nsites(tJDistributed const &block);
XDIAG_API int64_t dim(tJDistributed const &block);
XDIAG_API int64_t size(tJDistributed const &block);
XDIAG_API bool isreal(tJDistributed const &block);
XDIAG_API std::ostream &operator<<(std::ostream &out,
                                   tJDistributed const &block);
XDIAG_API std::string to_string(tJDistributed const &block);

class tJDistributedIterator {
public:
  tJDistributedIterator(tJDistributed const &block, bool begin);
  XDIAG_API tJDistributedIterator &operator++();
  XDIAG_API ProductState const &operator*() const;
  XDIAG_API bool operator!=(tJDistributedIterator const &rhs) const;

private:
  int64_t nsites_;
  mutable ProductState pstate_;
  basis::BasistJDistributedIterator it_;
};

} // namespace xdiag

#endif
