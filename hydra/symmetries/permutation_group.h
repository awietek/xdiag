#pragma once

#include <hydra/symmetries/permutation.h>
#include <hydra/symmetries/representation.h>
#include <vector>

namespace hydra {

class PermutationGroup {
public:
  PermutationGroup() = default;
  explicit PermutationGroup(std::vector<Permutation> const &permutations);

  inline int n_sites() const { return n_sites_; }
  inline int n_symmetries() const { return n_symmetries_; }
  inline int size() const { return n_symmetries_; }

  inline Permutation const &operator[](int sym) const {
    return permutations_[sym];
  }

  inline int inverse(int sym) const { return inverse_[sym]; }

  bool operator==(PermutationGroup const &rhs) const;
  bool operator!=(PermutationGroup const &rhs) const;

  PermutationGroup subgroup(std::vector<int> const &symmetry_numbers) const;

  inline auto begin() const { return permutations_.begin(); }
  inline auto end() const { return permutations_.end(); }

private:
  int n_sites_;
  int n_symmetries_;
  std::vector<Permutation> permutations_;
  std::vector<int> inverse_;
};

PermutationGroup allowed_subgroup(PermutationGroup const &group,
                                  Representation const &irrep);

} // namespace hydra
