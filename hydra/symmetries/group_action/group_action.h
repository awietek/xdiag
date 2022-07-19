#pragma once

#include <lila/external/gsl/span>

#include <hydra/symmetries/permutation_group.h>
#include <utility>

namespace hydra {

class GroupAction {
public:
  GroupAction() = default;
  GroupAction(PermutationGroup const &permutation_group);

  inline int n_sites() const { return n_sites_; }
  inline int n_symmetries() const { return n_symmetries_; }
  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }
  inline std::vector<int> const &permutation_array() const {
    return permutation_group_.permutation_array();
  }

  template <class bit_t> bit_t apply(int sym, bit_t state) const;
  template <class bit_t> bit_t representative(bit_t state) const;
  template <class bit_t>
  std::pair<bit_t, int> representative_index(bit_t state) const;
  template <class bit_t>
  std::pair<bit_t, gsl::span<int const>>
  representative_indices(bit_t state) const;
  template <class bit_t> double fermi_sign(int sym, bit_t state) const;
  template <class bit_t>
  std::vector<int> stabilizer_symmetries(bit_t state) const;

  bool operator==(GroupAction const &rhs) const;
  bool operator!=(GroupAction const &rhs) const;

private:
  int n_sites_;
  int n_symmetries_;
  PermutationGroup permutation_group_;

  mutable std::vector<int> indices_;
  mutable std::vector<int> fermi_work_;
};

} // namespace hydra
