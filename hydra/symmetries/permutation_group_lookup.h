#pragma once

#include <hydra/symmetries/permutation_group.h>
#include <tuple>

namespace hydra {

template <class bit_t> class PermutationGroupLookup {
public:
  PermutationGroupLookup() = default;
  PermutationGroupLookup(PermutationGroup const &permutation_group, int n_up);

  inline int n_sites() const { return n_sites_; }
  inline int n_symmetries() const { return n_symmetries_; }
  inline std::vector<int> const &permutation_array() const {
    return permutation_array_;
  }

  bit_t apply(int sym, bit_t state) const;
  bit_t representative(bit_t state) const;
  std::tuple<bit_t, int> representative_index(bit_t state) const;
  std::tuple<bit_t, int, const int *> representative_indices(bit_t state) const;
  double fermi_sign(int sym, bit_t state) const;
  std::vector<int> stabilizer_symmetries(bit_t state) const;

  bool operator==(PermutationGroupLookup const &rhs) const;
  bool operator!=(PermutationGroupLookup const &rhs) const;

private:
  int n_sites_;
  int n_up_;
  int n_symmetries_;
  std::vector<int> permutation_array_; // size = n_symmetries_*n_sites_

  std::vector<bit_t> translation_table_prefix_;
  std::vector<bit_t> translation_table_postfix_;

  mutable std::vector<int> indices_;
  mutable std::vector<int> fermi_work_;
};

} // namespace hydra
