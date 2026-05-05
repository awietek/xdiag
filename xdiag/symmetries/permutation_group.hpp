// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include <xdiag/armadillo.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Finite group of site permutations with precomputed inverse and multiplication
// table. Permutations are stored column-major in a (nsites x size) matrix so
// that each permutation's entries are contiguous in memory.
class PermutationGroup {
public:
  XDIAG_API PermutationGroup() = default;
  XDIAG_API explicit PermutationGroup(
      std::vector<Permutation> const &permutations);
  XDIAG_API explicit PermutationGroup(
      arma::Mat<int64_t> const &matrix); // nsites x n_permutations, cols=perms
  XDIAG_API PermutationGroup(int64_t *ptr, int64_t n_permutations,
                             int64_t nsites);

  XDIAG_API bool operator==(PermutationGroup const &rhs) const;
  XDIAG_API bool operator!=(PermutationGroup const &rhs) const;

  XDIAG_API int64_t size() const;
  XDIAG_API int64_t nsites() const;

  XDIAG_API Permutation operator[](int64_t sym) const;
  XDIAG_API int64_t const *ptr(int64_t sym) const; // pointer to sym-th column
  XDIAG_API int64_t inv(int64_t sym) const;
  XDIAG_API int64_t multiply(int64_t s1, int64_t s2) const;

  class iterator {
  public:
    iterator(PermutationGroup const *group, int64_t idx);
    Permutation operator*() const;
    iterator &operator++();
    bool operator!=(iterator const &rhs) const;

  private:
    PermutationGroup const *group_;
    int64_t idx_;
  };

  XDIAG_API iterator begin() const;
  XDIAG_API iterator end() const;

private:
  arma::Mat<int64_t> permutations_; // nsites x n_permutations
  std::vector<int64_t> inv_;
  arma::Mat<int64_t> multiply_;
};

XDIAG_API int64_t nsites(PermutationGroup const &group);
XDIAG_API int64_t size(PermutationGroup const &group);
XDIAG_API PermutationGroup subgroup(PermutationGroup const &group,
                                    std::vector<int64_t> const &symmetries);
XDIAG_API std::ostream &operator<<(std::ostream &out,
                                   PermutationGroup const &group);
XDIAG_API std::string to_string(PermutationGroup const &group);

} // namespace xdiag
