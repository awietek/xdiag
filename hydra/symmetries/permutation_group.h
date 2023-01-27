#pragma once

#include <hydra/symmetries/permutation.h>
#include <vector>

namespace hydra {

class PermutationGroup {
public:
  PermutationGroup() = default;
  explicit PermutationGroup(std::vector<Permutation> const &permutations);

  int n_sites() const;
  int n_symmetries() const;
  int size() const;
  Permutation const &operator[](int sym) const;
  int inverse(int sym) const;

  bool operator==(PermutationGroup const &rhs) const;
  bool operator!=(PermutationGroup const &rhs) const;

  PermutationGroup subgroup(std::vector<int> const &symmetry_numbers) const;

  using iterator_t = std::vector<Permutation>::const_iterator; 
  iterator_t begin() const;
  iterator_t end() const;

private:
  int n_sites_;
  int n_symmetries_;
  std::vector<Permutation> permutations_;
  std::vector<int> inverse_;
};

} // namespace hydra
