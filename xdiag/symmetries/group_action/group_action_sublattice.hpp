#pragma once
#include <array>
#include <utility>

#include <xdiag/bits/half_bit_types.hpp>
#include <xdiag/common.hpp>
#include <xdiag/extern/gsl/span>
#include <xdiag/symmetries/permutation_group.hpp>

namespace xdiag {

// Forward declaration for friend
namespace basis::spinhalf {
template <typename bit_t, int n_sublat> class BasisSublattice;
}

template <typename bit_t, int n_sublat> class GroupActionSublattice {
public:
  using half_bit_t = bits::half_bit_t<bit_t>;

  GroupActionSublattice() = default;
  GroupActionSublattice(PermutationGroup const &permutation_group);

  inline int64_t nsites() const { return nsites_; }
  inline int64_t n_symmetries() const { return n_symmetries_; }
  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }

  bit_t apply(int64_t sym, bit_t state) const;
  bit_t representative(bit_t state) const;
  std::pair<bit_t, int64_t> representative_sym(bit_t state) const;
  std::pair<bit_t, gsl::span<int64_t const>>
  representative_syms(bit_t state) const;

  bool operator==(GroupActionSublattice const &rhs) const;
  bool operator!=(GroupActionSublattice const &rhs) const;

// private:
  int64_t nsites_;
  int64_t n_symmetries_;
  PermutationGroup permutation_group_;

  int64_t nsites_sublat_;
  int64_t size_tables_;
  half_bit_t sublat_mask_;
  mutable std::vector<int64_t> representative_syms_;
  std::array<int, n_sublat> sublat_shift_;

  inline bit_t &sym_action(int sublat, int sym, half_bit_t bits) {
    return sym_action_[sublat][(int64_t)bits * n_symmetries_ + sym];
  }

  inline bit_t const &sym_action(int sublat, int sym, half_bit_t bits) const {
    return sym_action_[sublat][(int64_t)bits * n_symmetries_ + sym];
  }

  std::array<std::vector<half_bit_t>, n_sublat> reps_;
  std::array<std::vector<gsl::span<int64_t const>>, n_sublat> rep_syms_;
  std::array<std::vector<int64_t>, n_sublat> rep_syms_array_;
  std::array<std::vector<bit_t>, n_sublat> sym_action_;

  // friend class basis::spinhalf::BasisSublattice<bit_t, n_sublat>;
};

} // namespace xdiag
