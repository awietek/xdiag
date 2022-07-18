#pragma once

#include <lila/external/gsl/span>

#include <hydra/bitops/half_bit_types.h>
#include <hydra/symmetries/permutation_group.h>
#include <tuple>

namespace hydra {

template <typename bit_t, int n_sublat> class GroupActionSublattice {
public:
  using half_bit_t = bitops::half_bit_t<bit_t>;

  GroupActionSublattice() = default;
  GroupActionSublattice(PermutationGroup const &permutation_group);

  inline int n_sites() const { return n_sites_; }
  inline int n_symmetries() const { return n_symmetries_; }
  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }
  inline std::vector<int> const &permutation_array() const {
    return permutation_group_.permutation_array();
  }

  bit_t apply(int sym, bit_t state) const;
  bit_t representative(bit_t state) const;
  std::tuple<bit_t, int> representative_index(bit_t state) const;
  std::tuple<bit_t, int, const int *> representative_indices(bit_t state) const;

  bool operator==(GroupActionSublattice const &rhs) const;
  bool operator!=(GroupActionSublattice const &rhs) const;

private:
  int n_sites_;
  int n_symmetries_;
  PermutationGroup permutation_group_;

  int len_word_;
  int wordmask_;
  idx_t size_sublat_tables_;

  std::vector<std::pair<half_bit_t, gsl::span<int const>>> RepAndSyms[n_sublat];
  std::vector<int> RepSymsArray[n_sublat];
  std::vector<bit_t> SymAction[n_sublat];
};

} // namespace hydra
