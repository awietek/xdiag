#pragma once

#include <vector>

#include <hydra/symmetries/representation.h>
#include <lila/external/gsl/span>

namespace hydra {

class PermutationGroup {
public:
  PermutationGroup() = default;
  explicit PermutationGroup(std::vector<std::vector<int>> const &symmetries);
  PermutationGroup(int n_sites, int n_symmetries,
                   std::vector<int> const &permutation_array);
  PermutationGroup subgroup(std::vector<int> const &symmetry_numbers) const;

  inline int n_sites() const { return n_sites_; }
  inline int n_symmetries() const { return n_symmetries_; }
  inline std::vector<int> const &permutation_array() const {
    return permutation_array_;
  }

  inline int size() const { return n_symmetries_; }
  inline int permutation(int sym, int site) const {
    return permutation_array_[n_sites_ * sym + site];
  }
  inline gsl::span<int const> permutation(int sym) const{
    return gsl::span<int const>(permutation_array_.data() + n_sites_ * sym,
                                n_sites_);
  }

  bool operator==(PermutationGroup const &rhs) const;
  bool operator!=(PermutationGroup const &rhs) const;

private:
  int n_sites_;
  int n_symmetries_;
  std::vector<int> permutation_array_; // size = n_symmetries_*n_sites_
};

PermutationGroup allowed_subgroup(PermutationGroup const &group,
                                  Representation const &irrep);

} // namespace hydra
