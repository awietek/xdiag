#include "basis_symmetric_sz.hpp"

#include <xdiag/symmetries/operations/representative_list.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::basis::spinhalf {

template <class bit_t>
BasisSymmetricSz<bit_t>::BasisSymmetricSz(int64_t n_up,
                                          Representation const &irrep) try
    : n_sites_(irrep.group().n_sites()), n_up_(n_up),
      group_action_(irrep.group()), irrep_(irrep),
      combinations_indexing_(n_sites_, n_up) {
  check_n_sites_work_with_bits<bit_t>(n_sites_);

  if (n_up < 0) {
    XDIAG_THROW("Invalid value of nup: n_up < 0");
  } else if (n_up > n_sites_) {
    XDIAG_THROW("Invalid value of nup: n_up > n_sites");
  }

  if (isreal(irrep)) {
    arma::vec characters = irrep.characters().as<arma::vec>();
    std::tie(reps_, index_for_rep_, syms_, sym_limits_for_rep_, norms_) =
        symmetries::representatives_indices_symmetries_limits_norms<bit_t>(
            combinations_indexing_, group_action_, characters);
  } else {
    arma::cx_vec characters = irrep.characters().as<arma::cx_vec>();
    std::tie(reps_, index_for_rep_, syms_, sym_limits_for_rep_, norms_) =
        symmetries::representatives_indices_symmetries_limits_norms<bit_t>(
            combinations_indexing_, group_action_, characters);
  }

  size_ = (int64_t)reps_.size();
} catch (Error const &e) {
  XDIAG_RETHROW(e);
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

template class BasisSymmetricSz<uint32_t>;
template class BasisSymmetricSz<uint64_t>;

} // namespace xdiag::basis::spinhalf
