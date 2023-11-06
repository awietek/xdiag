#include "basis_symmetric_sz.h"

#include <hydra/symmetries/operations/representative_list.h>
#include <hydra/utils/logger.h>

namespace hydra::basis::spinhalf {

template <class bit_t>
BasisSymmetricSz<bit_t>::BasisSymmetricSz(int64_t n_sites, int64_t n_up,
                                          PermutationGroup group,
                                          Representation irrep)
    : n_sites_(n_sites), n_up_(n_up),
      group_action_(allowed_subgroup(group, irrep)), irrep_(irrep),
      combinations_indexing_(n_sites, n_up) {
  if ((n_up < 0) || (n_up > n_sites)) {
    throw(std::invalid_argument("Invalid value of nup"));
  } else if (n_sites < 0) {
    throw(std::invalid_argument("n_sites < 0"));
  } else if (n_sites != group.n_sites()) {
    throw(std::logic_error(
        "n_sites does not match the n_sites in PermutationGroup"));
  } else if (group_action_.n_symmetries() != irrep.size()) {
    throw(std::logic_error("PermutationGroup and Representation do not have "
                           "same number of elements"));
  }

  std::tie(reps_, index_for_rep_, syms_, sym_limits_for_rep_, norms_) =
      symmetries::representatives_indices_symmetries_limits_norms<bit_t>(
          combinations_indexing_, group_action_, irrep);

  size_ = (int64_t)reps_.size();
}

template <class bit_t>
typename std::vector<bit_t>::const_iterator
BasisSymmetricSz<bit_t>::begin() const {
  return reps_.begin();
}
template <class bit_t>
typename std::vector<bit_t>::const_iterator
BasisSymmetricSz<bit_t>::end() const {
  return reps_.end();
}

template <class bit_t> int64_t BasisSymmetricSz<bit_t>::dim() const {
  return size_;
}
template <class bit_t> int64_t BasisSymmetricSz<bit_t>::size() const {
  return size_;
}

template <class bit_t> int64_t BasisSymmetricSz<bit_t>::n_sites() const {
  return n_sites_;
}
template <class bit_t> int64_t BasisSymmetricSz<bit_t>::n_up() const {
  return n_up_;
}
template <class bit_t>
GroupActionLookup<bit_t> const &BasisSymmetricSz<bit_t>::group_action() const {
  return group_action_;
}
template <class bit_t>
Representation const &BasisSymmetricSz<bit_t>::irrep() const {
  return irrep_;
}

template <typename bit_t>
bool BasisSymmetricSz<bit_t>::operator==(
    BasisSymmetricSz<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_up_ == rhs.n_up_) &&
         (group_action_ == rhs.group_action_) && (irrep_ == rhs.irrep_);
}

template <typename bit_t>
bool BasisSymmetricSz<bit_t>::operator!=(
    BasisSymmetricSz<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class BasisSymmetricSz<uint16_t>;
template class BasisSymmetricSz<uint32_t>;
template class BasisSymmetricSz<uint64_t>;

} // namespace hydra::basis::spinhalf
