#include "basis_symmetric_no_sz.hpp"

#include <xdiag/symmetries/operations/representative_list.hpp>

namespace xdiag::basis::spinhalf {

template <class bit_t>
BasisSymmetricNoSz<bit_t>::BasisSymmetricNoSz(int64_t n_sites,
                                              PermutationGroup group,
                                              Representation irrep)
    : n_sites_(n_sites), group_action_(allowed_subgroup(group, irrep)),
      irrep_(irrep), subsets_basis_(n_sites) {

  if (n_sites < 0) {
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
          subsets_basis_, group_action_, irrep);
  size_ = (int64_t)reps_.size();
}

template <class bit_t>
typename std::vector<bit_t>::const_iterator
BasisSymmetricNoSz<bit_t>::begin() const {
  return reps_.begin();
}
  
template <class bit_t>
typename std::vector<bit_t>::const_iterator
BasisSymmetricNoSz<bit_t>::end() const {
  return reps_.end();
}
  
template <class bit_t> int64_t BasisSymmetricNoSz<bit_t>::dim() const {
  return size_;
}
  
template <class bit_t> int64_t BasisSymmetricNoSz<bit_t>::size() const {
  return size_;
}

template <class bit_t> int64_t BasisSymmetricNoSz<bit_t>::n_sites() const {
  return n_sites_;
}
template <class bit_t>
GroupActionLookup<bit_t> const &
BasisSymmetricNoSz<bit_t>::group_action() const {
  return group_action_;
}
template <class bit_t>
Representation const &BasisSymmetricNoSz<bit_t>::irrep() const {
  return irrep_;
}

template <typename bit_t>
bool BasisSymmetricNoSz<bit_t>::operator==(
    BasisSymmetricNoSz<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (group_action_ == rhs.group_action_) &&
         (irrep_ == rhs.irrep_);
}

template <typename bit_t>
bool BasisSymmetricNoSz<bit_t>::operator!=(
    BasisSymmetricNoSz<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class BasisSymmetricNoSz<uint32_t>;
template class BasisSymmetricNoSz<uint64_t>;

} // namespace xdiag::basis::spinhalf
