#include "basis_symmetric_np.hpp"

#include <xdiag/combinatorics/combinations_indexing.hpp>
#include <xdiag/combinatorics/fermi_table.hpp>

#include <xdiag/symmetries/operations/group_action_operations.hpp>
#include <xdiag/symmetries/operations/representative_list.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>

namespace xdiag::basis::electron {

template <class bit_t>
BasisSymmetricNp<bit_t>::BasisSymmetricNp(int64_t n_sites, int64_t nup,
                                          int64_t ndn, PermutationGroup group,
                                          Representation irrep)
    : n_sites_(n_sites), n_up_(nup), n_dn_(ndn),
      group_action_(allowed_subgroup(group, irrep)), irrep_(irrep),
      raw_ups_size_(combinatorics::binomial(n_sites, nup)),
      raw_dns_size_(combinatorics::binomial(n_sites, ndn)),
      lintable_ups_(n_sites, nup), lintable_dns_(n_sites, ndn),
      fermi_table_ups_(n_sites, nup, allowed_subgroup(group, irrep)),
      fermi_table_dns_(n_sites, ndn, allowed_subgroup(group, irrep)) {

  using combinatorics::Combinations;
  if (n_sites < 0) {
    throw(std::invalid_argument("n_sites < 0"));
  } else if ((nup < 0) || (nup > n_sites)) {
    throw(std::invalid_argument("Invalid value of nup"));
  } else if ((ndn < 0) || (ndn > n_sites)) {
    throw(std::invalid_argument("Invalid value of ndn"));
  } else if (n_sites != group.n_sites()) {
    throw(std::logic_error(
        "n_sites does not match the n_sites in PermutationGroup"));
  }

  std::tie(reps_up_, idces_up_, syms_up_, sym_limits_up_) =
      symmetries::representatives_indices_symmetries_limits<bit_t>(
          combinatorics::CombinationsIndexing<bit_t>(n_sites, nup),
          group_action_);
  std::tie(dns_storage_, norms_storage_, dns_limits_, ups_offset_, size_) =
      symmetries::electron_dns_norms_limits_offset_size(
          reps_up_, Combinations<bit_t>(n_sites, ndn), group_action_, irrep_);
}

template <class bit_t> int64_t BasisSymmetricNp<bit_t>::n_sites() const {
  return n_sites_;
}
template <class bit_t> int64_t BasisSymmetricNp<bit_t>::n_up() const {
  return n_up_;
}
template <class bit_t> int64_t BasisSymmetricNp<bit_t>::n_dn() const {
  return n_dn_;
}

template <class bit_t> int64_t BasisSymmetricNp<bit_t>::dim() const {
  return size_;
}
template <class bit_t> int64_t BasisSymmetricNp<bit_t>::size() const {
  return size_;
}
template <typename bit_t>
typename BasisSymmetricNp<bit_t>::iterator_t
BasisSymmetricNp<bit_t>::begin() const {
  return iterator_t(*this, true);
}
template <typename bit_t>
typename BasisSymmetricNp<bit_t>::iterator_t
BasisSymmetricNp<bit_t>::end() const {
  return iterator_t(*this, false);
}

template <class bit_t>
GroupActionLookup<bit_t> const &BasisSymmetricNp<bit_t>::group_action() const {
  return group_action_;
}
template <class bit_t>
Representation const &BasisSymmetricNp<bit_t>::irrep() const {
  return irrep_;
}

template class BasisSymmetricNp<uint32_t>;
template class BasisSymmetricNp<uint64_t>;

template <typename bit_t>
BasisSymmetricNpIterator<bit_t>::BasisSymmetricNpIterator(
    BasisSymmetricNp<bit_t> const &basis, bool begin)
    : basis_(basis), up_idx_(begin ? 0 : basis.n_rep_ups()), dn_idx_(0) {
  if (basis.dim() == 0) {
    up_idx_ = 0;
  }

  if ((basis.n_rep_ups() > 0) && begin) {
    dns_for_ups_rep_ = basis.dns_for_ups_rep(basis.rep_ups(0));
  }
}

template <typename bit_t>
BasisSymmetricNpIterator<bit_t> &BasisSymmetricNpIterator<bit_t>::operator++() {
  ++dn_idx_;
  if (dn_idx_ == dns_for_ups_rep_.size()) {
    dn_idx_ = 0;
    do {
      ++up_idx_;
      if (up_idx_ == basis_.n_rep_ups()) {
        return *this;
      }
      bit_t ups = basis_.rep_ups(up_idx_);
      dns_for_ups_rep_ = basis_.dns_for_ups_rep(ups);
    } while (dns_for_ups_rep_.size() == 0);
  }
  return *this;
}

template <typename bit_t>
std::pair<bit_t, bit_t> BasisSymmetricNpIterator<bit_t>::operator*() const {
  return {basis_.rep_ups(up_idx_), dns_for_ups_rep_[dn_idx_]};
}

template <typename bit_t>
bool BasisSymmetricNpIterator<bit_t>::operator!=(
    BasisSymmetricNpIterator<bit_t> const &rhs) const {
  return (up_idx_ != rhs.up_idx_) || (dn_idx_ != rhs.dn_idx_);
}

template class BasisSymmetricNpIterator<uint32_t>;
template class BasisSymmetricNpIterator<uint64_t>;

} // namespace xdiag::basis::electron
