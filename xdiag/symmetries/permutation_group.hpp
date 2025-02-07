#pragma once

#include <ostream>
#include <string>
#include <vector>

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/symmetries/permutation.hpp>

namespace xdiag {

class PermutationGroup {
public:
  using iterator_t = std::vector<Permutation>::const_iterator;

  XDIAG_API PermutationGroup() = default;
  XDIAG_API explicit PermutationGroup(
      std::vector<Permutation> const &permutations);
  XDIAG_API explicit PermutationGroup(arma::Mat<int64_t> const &matrix);
  XDIAG_API PermutationGroup(int64_t *ptr, int64_t n_permutations,
                             int64_t nsites);

  XDIAG_API bool operator==(PermutationGroup const &rhs) const;
  XDIAG_API bool operator!=(PermutationGroup const &rhs) const;

  XDIAG_API int64_t size() const;
  XDIAG_API int64_t nsites() const;

  XDIAG_API Permutation const &operator[](int64_t sym) const;
  XDIAG_API int64_t inv(int64_t sym) const;
  XDIAG_API int64_t multiply(int64_t s1, int64_t s2) const;

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;

private:
  std::vector<Permutation> permutations_;
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
