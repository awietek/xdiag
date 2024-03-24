#pragma once

#include <hydra/extern/gsl/span>

#include <hydra/symmetries/permutation_group.h>
#include <utility>

namespace hydra {

class GroupAction {
public:
  GroupAction() = default;
  GroupAction(PermutationGroup const &permutation_group);

  inline int64_t n_sites() const { return n_sites_; }
  inline int64_t n_symmetries() const { return n_symmetries_; }
  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }

  template <class bit_t> bit_t apply(int64_t sym, bit_t state) const;
  template <class bit_t> bit_t representative(bit_t state) const;
  template <class bit_t>
  std::pair<bit_t, int64_t> representative_sym(bit_t state) const;
  template <class bit_t>
  std::pair<bit_t, gsl::span<int64_t const>> representative_syms(bit_t state) const;

  template <class bit_t> double fermi_sign(int64_t sym, bit_t state) const;
  template <class bit_t>
  std::vector<int64_t> stabilizer_symmetries(bit_t state) const;

  bool operator==(GroupAction const &rhs) const;
  bool operator!=(GroupAction const &rhs) const;

private:
  int64_t n_sites_;
  int64_t n_symmetries_;
  PermutationGroup permutation_group_;

  mutable std::vector<int64_t> indices_;
  mutable std::vector<int64_t> fermi_work_;
};

} // namespace hydra
