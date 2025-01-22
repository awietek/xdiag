#include "group_action.hpp"

#include <xdiag/symmetries/operations/fermi_sign.hpp>
#include <xdiag/symmetries/operations/group_action_operations.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>

namespace xdiag {

GroupAction::GroupAction(PermutationGroup const &permutation_group)
    : nsites_(permutation_group.nsites()),
      n_symmetries_(permutation_group.size()),
      permutation_group_(permutation_group), indices_(n_symmetries_, 0),
      fermi_work_(2 * permutation_group.nsites(), 0) {}

template <class bit_t> bit_t GroupAction::apply(int64_t sym, bit_t state) const {
  return permutation_group_[sym].apply(state);
}
template uint16_t GroupAction::apply<uint16_t>(int64_t, uint16_t) const;
template uint32_t GroupAction::apply<uint32_t>(int64_t, uint32_t) const;
template uint64_t GroupAction::apply<uint64_t>(int64_t, uint64_t) const;

template <class bit_t> bit_t GroupAction::representative(bit_t state) const {
  bit_t rep = std::numeric_limits<bit_t>::max();

  for (auto const &p : permutation_group_) {
    bit_t trans = p.apply(state);
    if (trans < rep) {
      rep = trans;
    }
  }
  return rep;
}
template uint16_t GroupAction::representative<uint16_t>(uint16_t) const;
template uint32_t GroupAction::representative<uint32_t>(uint32_t) const;
template uint64_t GroupAction::representative<uint64_t>(uint64_t) const;

template <class bit_t>
std::pair<bit_t, int64_t> GroupAction::representative_sym(bit_t state) const {
  bit_t rep = std::numeric_limits<bit_t>::max();
  int64_t idx = 0;

  int64_t sym=0;
  for (auto const &p : permutation_group_) {
    bit_t trans = p.apply(state);
    if (trans < rep) {
      rep = trans;
      idx = sym;
    }
    ++sym;
  }
  return {rep, idx};
}
template std::pair<uint16_t, int64_t>
    GroupAction::representative_sym<uint16_t>(uint16_t) const;
template std::pair<uint32_t, int64_t>
    GroupAction::representative_sym<uint32_t>(uint32_t) const;
template std::pair<uint64_t, int64_t>
    GroupAction::representative_sym<uint64_t>(uint64_t) const;

template <class bit_t>
std::pair<bit_t, gsl::span<int64_t const>>
GroupAction::representative_syms(bit_t state) const {
  bit_t rep = std::numeric_limits<bit_t>::max();
  gsl::span<int64_t const>::size_type n_indices = 0;

  int64_t sym=0;
  for (auto const &p : permutation_group_) {
    bit_t trans = p.apply(state);
    if (trans < rep) {
      rep = trans;
      n_indices = 1;
      indices_[0] = sym;
    } else if (trans == rep) {
      indices_[n_indices++] = sym;
    }
    ++sym;
  }
  return {rep, {indices_.data(), n_indices}};
}

template std::pair<uint16_t, gsl::span<int64_t const>>
    GroupAction::representative_syms<uint16_t>(uint16_t) const;
template std::pair<uint32_t, gsl::span<int64_t const>>
    GroupAction::representative_syms<uint32_t>(uint32_t) const;
template std::pair<uint64_t, gsl::span<int64_t const>>
    GroupAction::representative_syms<uint64_t>(uint64_t) const;

template <class bit_t>
double GroupAction::fermi_sign(int64_t sym, bit_t state) const {
  return symmetries::fermi_sign_of_permutation(state, permutation_group_[sym],
                                               fermi_work_);
}
template double GroupAction::fermi_sign<uint16_t>(int64_t, uint16_t) const;
template double GroupAction::fermi_sign<uint32_t>(int64_t, uint32_t) const;
template double GroupAction::fermi_sign<uint64_t>(int64_t, uint64_t) const;

template <class bit_t>
std::vector<int64_t> GroupAction::stabilizer_symmetries(bit_t bits) const {
  return symmetries::stabilizer_symmetries(bits, *this);
}
template std::vector<int64_t>
    GroupAction::stabilizer_symmetries<uint16_t>(uint16_t) const;
template std::vector<int64_t>
    GroupAction::stabilizer_symmetries<uint32_t>(uint32_t) const;
template std::vector<int64_t>
    GroupAction::stabilizer_symmetries<uint64_t>(uint64_t) const;

bool GroupAction::operator==(GroupAction const &rhs) const {
  return (nsites_ == rhs.nsites_) && (n_symmetries_ == rhs.n_symmetries_) &&
         (permutation_group_ == rhs.permutation_group_);
}

bool GroupAction::operator!=(GroupAction const &rhs) const {
  return !operator==(rhs);
}

} // namespace xdiag
