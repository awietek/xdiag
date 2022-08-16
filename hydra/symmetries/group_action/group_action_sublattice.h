#pragma once

#include <lila/external/gsl/span>

#include <array>
#include <hydra/bitops/half_bit_types.h>
#include <hydra/symmetries/permutation_group.h>
#include <utility>

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

  bit_t apply(int sym, bit_t state) const;
  bit_t representative(bit_t state) const;
  std::pair<bit_t, int> representative_sym(bit_t state) const;
  std::pair<bit_t, gsl::span<int const>> representative_syms(bit_t state) const;

  bool operator==(GroupActionSublattice const &rhs) const;
  bool operator!=(GroupActionSublattice const &rhs) const;

private:
  int n_sites_;
  int n_symmetries_;
  PermutationGroup permutation_group_;

  int n_sites_sublat_;
  idx_t size_tables_;
  half_bit_t sublat_mask_;
  mutable std::vector<int> representative_syms_;
  std::array<int, n_sublat> sublat_shift_;

  inline bit_t &sym_action(int sublat, int sym, half_bit_t bits) {
    return sym_action_[sublat][(idx_t)bits * n_symmetries_ + sym];
  }

  inline bit_t const &sym_action(int sublat, int sym, half_bit_t bits) const {
    return sym_action_[sublat][(idx_t)bits * n_symmetries_ + sym];
  }

  std::array<std::vector<half_bit_t>, n_sublat> reps_;
  std::array<std::vector<gsl::span<int const>>, n_sublat> rep_syms_;
  std::array<std::vector<int>, n_sublat> rep_syms_array_;
  std::array<std::vector<bit_t>, n_sublat> sym_action_;
};

} // namespace hydra
