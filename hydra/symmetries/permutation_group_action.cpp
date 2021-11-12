#include "permutation_group_action.h"

#include <hydra/symmetries/fermi_sign.h>
#include <hydra/symmetries/symmetry_operations.h>

namespace hydra {

PermutationGroupAction::PermutationGroupAction(
    PermutationGroup const &permutation_group)
    : n_sites_(permutation_group.n_sites()),
      n_symmetries_(permutation_group.n_symmetries()),
      permutation_group_(permutation_group), indices_(n_symmetries_, 0),
      fermi_work_(2 * permutation_group.n_sites(), 0) {}

template <class bit_t>
bit_t PermutationGroupAction::apply(int sym, bit_t state) const {
  return symmetries::apply_permutation(
      state, n_sites_, permutation_array().data() + sym * n_sites_);
}
template uint16_t PermutationGroupAction::apply<uint16_t>(int, uint16_t) const;
template uint32_t PermutationGroupAction::apply<uint32_t>(int, uint32_t) const;
template uint64_t PermutationGroupAction::apply<uint64_t>(int, uint64_t) const;

template <class bit_t>
bit_t PermutationGroupAction::representative(bit_t state) const {
  bit_t rep = std::numeric_limits<bit_t>::max();
  const int *sym_ptr = permutation_array().data();
  for (int sym = 0; sym < n_symmetries_; ++sym) {
    bit_t trans = symmetries::apply_permutation(state, n_sites_, sym_ptr);
    if (trans < rep)
      rep = trans;
    sym_ptr += n_sites_;
  }
  return rep;
}
template uint16_t
    PermutationGroupAction::representative<uint16_t>(uint16_t) const;
template uint32_t
    PermutationGroupAction::representative<uint32_t>(uint32_t) const;
template uint64_t
    PermutationGroupAction::representative<uint64_t>(uint64_t) const;

template <class bit_t>
std::tuple<bit_t, int>
PermutationGroupAction::representative_index(bit_t state) const {
  bit_t rep = std::numeric_limits<bit_t>::max();
  int idx = 0;
  const int *sym_ptr = permutation_array().data();
  for (int sym = 0; sym < n_symmetries_; ++sym) {
    bit_t trans = symmetries::apply_permutation(state, n_sites_, sym_ptr);
    if (trans < rep) {
      rep = trans;
      idx = sym;
    }
    sym_ptr += n_sites_;
  }
  return {rep, idx};
}
template std::tuple<uint16_t, int>
    PermutationGroupAction::representative_index<uint16_t>(uint16_t) const;
template std::tuple<uint32_t, int>
    PermutationGroupAction::representative_index<uint32_t>(uint32_t) const;
template std::tuple<uint64_t, int>
    PermutationGroupAction::representative_index<uint64_t>(uint64_t) const;

template <class bit_t>
std::tuple<bit_t, int, const int *>
PermutationGroupAction::representative_indices(bit_t state) const {
  bit_t rep = std::numeric_limits<bit_t>::max();
  int n_indices = 0;
  const int *sym_ptr = permutation_array().data();
  for (int sym = 0; sym < n_symmetries_; ++sym) {
    bit_t trans = symmetries::apply_permutation(state, n_sites_, sym_ptr);
    if (trans < rep) {
      rep = trans;
      n_indices = 1;
      indices_[0] = sym;
    } else if (trans == rep) {
      indices_[n_indices++] = sym;
    }
    sym_ptr += n_sites_;
  }
  return {rep, n_indices, indices_.data()};
}

template std::tuple<uint16_t, int, const int *>
    PermutationGroupAction::representative_indices<uint16_t>(uint16_t) const;
template std::tuple<uint32_t, int, const int *>
    PermutationGroupAction::representative_indices<uint32_t>(uint32_t) const;
template std::tuple<uint64_t, int, const int *>
    PermutationGroupAction::representative_indices<uint64_t>(uint64_t) const;

template <class bit_t>
double PermutationGroupAction::fermi_sign(int sym, bit_t state) const {
  return symmetries::fermi_sign_of_permutation(
      state, permutation_array().data() + sym * n_sites_, fermi_work_.data());
}
template double PermutationGroupAction::fermi_sign<uint16_t>(int,
                                                             uint16_t) const;
template double PermutationGroupAction::fermi_sign<uint32_t>(int,
                                                             uint32_t) const;
template double PermutationGroupAction::fermi_sign<uint64_t>(int,
                                                             uint64_t) const;

template <class bit_t>
std::vector<int>
PermutationGroupAction::stabilizer_symmetries(bit_t bits) const {
  return symmetries::stabilizer_symmetries(bits, *this);
}
template std::vector<int>
    PermutationGroupAction::stabilizer_symmetries<uint16_t>(uint16_t) const;
template std::vector<int>
    PermutationGroupAction::stabilizer_symmetries<uint32_t>(uint32_t) const;
template std::vector<int>
    PermutationGroupAction::stabilizer_symmetries<uint64_t>(uint64_t) const;

bool PermutationGroupAction::operator==(
    PermutationGroupAction const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_symmetries_ == rhs.n_symmetries_) &&
         (permutation_group_ == rhs.permutation_group_);
}

bool PermutationGroupAction::operator!=(
    PermutationGroupAction const &rhs) const {
  return !operator==(rhs);
}

} // namespace hydra
