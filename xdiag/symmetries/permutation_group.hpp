#pragma once
#include <vector>

#include <xdiag/io/toml/file_toml_handler.hpp>
#include <xdiag/symmetries/permutation.hpp>

namespace xdiag {

class PermutationGroup {
public:
  PermutationGroup() = default;
  explicit PermutationGroup(std::vector<Permutation> const &permutations);

  int64_t n_sites() const;
  int64_t n_symmetries() const;
  int64_t size() const;
  Permutation const &operator[](int64_t sym) const;
  int64_t inverse(int64_t sym) const;

  bool operator==(PermutationGroup const &rhs) const;
  bool operator!=(PermutationGroup const &rhs) const;
  operator bool() const;

  PermutationGroup subgroup(std::vector<int64_t> const &symmetry_numbers) const;
  using iterator_t = std::vector<Permutation>::const_iterator;
  iterator_t begin() const;
  iterator_t end() const;

private:
  int64_t n_sites_;
  int64_t n_symmetries_;
  std::vector<Permutation> permutations_;
  std::vector<int64_t> inverse_;
};

std::ostream &operator<<(std::ostream &out, PermutationGroup const &group);
std::string to_string(PermutationGroup const &group);

} // namespace xdiag
