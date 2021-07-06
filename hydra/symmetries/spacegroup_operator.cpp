#include "spacegroup_operator.h"

#include <lila/utils/logger.h>
#include <limits>

namespace hydra {

template <class bit_t>
SpaceGroupOperator<bit_t>::SpaceGroupOperator(
    int n_sites, std::vector<int> const &permutation_array)
    : n_sites_(n_sites), permutation_array_(permutation_array),
      n_sym_(permutation_array.size() / n_sites), indices_(n_sym_, 0),
      fermi_work_(2 * n_sites, 0) {
  if (n_sites <= 0)
    lila::Log.err("Error constructing SpaceGroupOperator: "
                  "invalid n_sites");

  ;

  if (permutation_array.size() % n_sites != 0)
    lila::Log.err("Error constructing SpaceGroupOperator: "
                  "size of permutation_array not a multiple of n_sites");
}

template <class bit_t>
bit_t SpaceGroupOperator<bit_t>::representative(bit_t state) const {
  bit_t rep = std::numeric_limits<bit_t>::max();
  const int *sym_ptr = permutation_array_.data();
  for (int sym = 0; sym < n_sym_; ++sym) {
    bit_t trans = symmetries::apply_permutation(state, n_sites_, sym_ptr);
    if (trans < rep)
      rep = trans;
    sym_ptr += n_sites_;
  }
  return rep;
}

template <class bit_t>
std::tuple<bit_t, int>
SpaceGroupOperator<bit_t>::representative_index(bit_t state) const {
  bit_t rep = std::numeric_limits<bit_t>::max();
  int idx = 0;
  const int *sym_ptr = permutation_array_.data();
  for (int sym = 0; sym < n_sym_; ++sym) {
    bit_t trans = symmetries::apply_permutation(state, n_sites_, sym_ptr);
    if (trans < rep) {
      rep = trans;
      idx = sym;
    }
    sym_ptr += n_sites_;
  }
  return {rep, idx};
}

template <class bit_t>
std::tuple<bit_t, int, const int *>
SpaceGroupOperator<bit_t>::representative_indices(bit_t state) const {
  bit_t rep = std::numeric_limits<bit_t>::max();
  int n_indices = 0;
  const int *sym_ptr = permutation_array_.data();
  for (int sym = 0; sym < n_sym_; ++sym) {
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

template <class bit_t>
bool SpaceGroupOperator<bit_t>::operator==(
    SpaceGroupOperator const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_sym_ == rhs.n_sym_) &&
         (permutation_array_ == rhs.permutation_array_);
}

template <class bit_t>
bool SpaceGroupOperator<bit_t>::operator!=(
    SpaceGroupOperator const &rhs) const {
  return !operator==(rhs);
}

template class SpaceGroupOperator<uint16>;
template class SpaceGroupOperator<uint32>;
template class SpaceGroupOperator<uint64>;

} // namespace hydra
