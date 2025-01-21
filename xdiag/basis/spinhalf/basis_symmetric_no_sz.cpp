#include "basis_symmetric_no_sz.hpp"

#include <xdiag/symmetries/operations/representative_list.hpp>

namespace xdiag::basis::spinhalf {

template <class bit_t>
BasisSymmetricNoSz<bit_t>::BasisSymmetricNoSz(Representation const &irrep) try
    : n_sites_(irrep.group().n_sites()), group_action_(irrep.group()),
      irrep_(irrep), subsets_basis_(n_sites_) {
  check_n_sites_work_with_bits<bit_t>(n_sites_);

  if (isreal(irrep)) {
    arma::vec characters = irrep.characters().as<arma::vec>();
    std::tie(reps_, index_for_rep_, syms_, sym_limits_for_rep_, norms_) =
        symmetries::representatives_indices_symmetries_limits_norms<bit_t>(
            subsets_basis_, group_action_, characters);
  } else {
    arma::cx_vec characters = irrep.characters().as<arma::cx_vec>();
    std::tie(reps_, index_for_rep_, syms_, sym_limits_for_rep_, norms_) =
        symmetries::representatives_indices_symmetries_limits_norms<bit_t>(
            subsets_basis_, group_action_, characters);
  }

  size_ = (int64_t)reps_.size();
} catch (Error const &e) {
  XDIAG_RETHROW(e);
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
