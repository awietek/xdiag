// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "permutation_group.hpp"

#include <algorithm>
#include <xdiag/common.hpp>

namespace xdiag {

static void
check_valid_group(std::vector<Permutation> const &permutations) try {
  if (permutations.size() == 0) {
    XDIAG_THROW(
        "A Permutation group must consist of at least one Permutation.");
  }

  int64_t nsites = permutations[0].size();

  // Check whether all permutations have same number of sites
  for (auto p : permutations) {
    if (p.size() != nsites) {
      XDIAG_THROW("Not all Permutations have the same number of sites");
    }
  }

  // Check whether every permutation is unique
  auto a = permutations;
  std::sort(a.begin(), a.end(),
            [](Permutation const &p1, Permutation const &p2) {
              auto const &a1 = p1.array();
              auto const &a2 = p2.array();
              return std::lexicographical_compare(a1.begin(), a1.end(),
                                                  a2.begin(), a2.end());
            });
  auto it = std::adjacent_find(a.begin(), a.end());
  if (it != a.end()) {
    XDIAG_THROW(fmt::format("Non-unique Permutation found. The following "
                            "permutation seems to appear twice or more:\n{}",
                            to_string(*it)));
  }

  // Check whether identity is contained
  auto id = Permutation(nsites);
  if (std::find(permutations.begin(), permutations.end(), id) ==
      permutations.end()) {
    XDIAG_THROW("Identity element not found");
  }
}
XDIAG_CATCH

static std::vector<int64_t>
compute_inverse(std::vector<Permutation> const &permutations) try {
  std::vector<int64_t> inv(permutations.size());
  int64_t idx = 0;
  for (auto p : permutations) {
    auto pinv = p.inv();
    auto it = std::find(permutations.begin(), permutations.end(), pinv);
    if (it == permutations.end()) {
      XDIAG_THROW(fmt::format("Inverse element not found. The inverse of the "
                              "following Permutation has not been found\n{}",
                              to_string(p)));
    } else {
      inv[idx] = std::distance(permutations.begin(), it);
    }
    idx++;
  }
  return inv;
}
XDIAG_CATCH

static arma::Mat<int64_t>
compute_multiply(std::vector<Permutation> const &permutations) try {
  int64_t n = permutations.size();
  arma::Mat<int64_t> multiply(n, n);
  int64_t i1 = 0;
  for (auto p1 : permutations) {
    int64_t i2 = 0;
    for (auto p2 : permutations) {
      auto p = p1 * p2;
      auto it = std::find(permutations.begin(), permutations.end(), p);
      if (it == permutations.end()) {
        XDIAG_THROW(fmt::format(
            "Group multiplication not closed. The product of the following "
            "two permutations is not an element of the group:\n{}\n{}",
            to_string(p1), to_string(p2)));
      } else {
        multiply(i1, i2) = std::distance(permutations.begin(), it);
      }
      ++i2;
    }
    ++i1;
  }
  return multiply;
}
XDIAG_CATCH

PermutationGroup::PermutationGroup(
    std::vector<Permutation> const &permutations) try
    : permutations_(permutations), inv_(compute_inverse(permutations)),
      multiply_(compute_multiply(permutations)) {
  check_valid_group(permutations);
}
XDIAG_CATCH

static std::vector<Permutation>
permutations_from_matrix(arma::Mat<int64_t> const &matrix) try {
  std::vector<Permutation> permutations;
  for (int64_t i = 0; i < matrix.n_rows; ++i) {
    auto perm = Permutation(arma::Col<int64_t>(matrix.row(i).t()));
    permutations.push_back(perm);
  }
  return permutations;
}
XDIAG_CATCH

PermutationGroup::PermutationGroup(arma::Mat<int64_t> const &matrix) try
    : PermutationGroup(permutations_from_matrix(matrix)) {}
XDIAG_CATCH

PermutationGroup::PermutationGroup(int64_t *ptr, int64_t n_permutations,
                                   int64_t nsites) try
    : PermutationGroup(permutations_from_matrix(
          arma::Mat<int64_t>(ptr, n_permutations, nsites, true))) {}
XDIAG_CATCH

int64_t PermutationGroup::nsites() const {
  return (permutations_.size() == 0) ? 0 : permutations_[0].size();
}
int64_t PermutationGroup::size() const { return permutations_.size(); }
Permutation const &PermutationGroup::operator[](int64_t sym) const try {
  if ((sym < 0) || (sym >= size())) {
    XDIAG_THROW("Invalid symmetry index");
  } else {
    return permutations_[sym];
  }
}
XDIAG_CATCH

int64_t PermutationGroup::inv(int64_t sym) const try {
  if ((sym < 0) || (sym >= size())) {
    XDIAG_THROW("Invalid symmetry index");
  } else {
    return inv_[sym];
  }
}
XDIAG_CATCH

int64_t PermutationGroup::multiply(int64_t s1, int64_t s2) const try {
  if ((s1 < 0) || (s1 >= size()) || (s2 < 0) || (s2 >= size())) {
    XDIAG_THROW("Invalid symmetry index");
  } else {
    return multiply_(s1, s2);
  }
}
XDIAG_CATCH

int64_t nsites(PermutationGroup const &group) { return group.nsites(); }
int64_t size(PermutationGroup const &group) { return group.size(); }
PermutationGroup subgroup(PermutationGroup const &group,
                          std::vector<int64_t> const &symmetries) try {
  std::vector<Permutation> subgroup_permutations;
  for (int64_t n_sym : symmetries) {

    if ((0 > n_sym) || (n_sym >= group.size())) {
      XDIAG_THROW("Invalid symmetry index.");
    }
    subgroup_permutations.push_back(group[n_sym]);
  }
  return PermutationGroup(subgroup_permutations);
}
XDIAG_CATCH

bool PermutationGroup::operator==(PermutationGroup const &rhs) const {
  return (permutations_ == rhs.permutations_);
}

bool PermutationGroup::operator!=(PermutationGroup const &rhs) const {
  return !operator==(rhs);
}

PermutationGroup::iterator_t PermutationGroup::begin() const {
  return permutations_.begin();
}
PermutationGroup::iterator_t PermutationGroup::end() const {
  return permutations_.end();
}

std::ostream &operator<<(std::ostream &out, PermutationGroup const &group) {
  out << "nsites  : " << group.nsites() << "\n";
  out << "size     : " << group.size() << "\n";
  for (auto const &p : group) {
    out << p << "\n";
  }
  return out;
}
std::string to_string(PermutationGroup const &group) {
  return to_string_generic(group);
}

} // namespace xdiag
