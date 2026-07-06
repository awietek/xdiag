// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include <xdiag/armadillo.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Permutation of n sites: bijection {0,...,n-1} -> {0,...,n-1} stored as p[i].
// Composition: (p1 * p2)[i] = p1[p2[i]], i.e. p2 is applied first.
class XDIAG_API Permutation {
public:
  Permutation() = default;
  explicit Permutation(int64_t size); // creates identity permutation
  Permutation(std::initializer_list<int64_t> list);
  explicit Permutation(std::vector<int32_t> const &array);
  explicit Permutation(std::vector<int64_t> const &array);
  explicit Permutation(arma::Col<int64_t> const &array);
  Permutation(int64_t const *ptr, int64_t size);

  int64_t size() const;
  int64_t operator[](int64_t i) const;
  Permutation inv() const;
  std::vector<int64_t> const &array() const;

  Permutation &operator*=(Permutation const &rhs);
  bool operator==(Permutation const &rhs) const;
  bool operator!=(Permutation const &rhs) const;

private:
  std::vector<int64_t> array_;
};

XDIAG_API int64_t size(Permutation const &p);
XDIAG_API Permutation multiply(Permutation const &p1, Permutation const &p2);
XDIAG_API Permutation operator*(Permutation const &p1, Permutation const &p2);
XDIAG_API Permutation inv(Permutation const &p);
XDIAG_API Permutation pow(Permutation const &p, int64_t power);

XDIAG_API std::ostream &operator<<(std::ostream &out, Permutation const &perm);
XDIAG_API std::string to_string(Permutation const &perm);

} // namespace xdiag
