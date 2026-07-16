// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>
#include <xdiag/armadillo.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// A ProductState is a list of per-site local quantum-number indices. Each entry
// is an integer in [0, d_i) where d_i is the local Hilbert space dimension at
// site i (e.g. spinhalf: Dn=0, Up=1; boson: the occupation number). The mapping
// from these integers to human-readable labels is a display concern handled by
// the block layer (see to_string(ProductState, Block)).
class XDIAG_API ProductState {
public:
  using iterator_t = std::vector<int64_t>::const_iterator;

  ProductState() = default;
  explicit ProductState(int64_t nsites);
  ProductState(std::initializer_list<int64_t> list);
  explicit ProductState(std::vector<int32_t> const &local_states);
  explicit ProductState(std::vector<int64_t> const &local_states);
  explicit ProductState(arma::Col<int64_t> const &local_states);

  int64_t operator[](int64_t i) const;
  int64_t &operator[](int64_t i);

  int64_t size() const;
  int64_t nsites() const;
  void push_back(int64_t l);

  iterator_t begin() const;
  iterator_t end() const;

  bool operator==(ProductState const &rhs) const;
  bool operator!=(ProductState const &rhs) const;

private:
  std::vector<int64_t> local_states_;
};

XDIAG_API int64_t size(ProductState const &p);
XDIAG_API int64_t nsites(ProductState const &p);
XDIAG_API std::ostream &operator<<(std::ostream &out,
                                   ProductState const &state);
XDIAG_API std::string to_string(ProductState const &state);

} // namespace xdiag
