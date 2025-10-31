// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_symmetric_no_np.hpp"

#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/symmetries/operations/group_action_operations.hpp>
#include <xdiag/symmetries/operations/representative_list.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>

namespace xdiag::basis::electron {

template <class bit_t>
BasisSymmetricNoNp<bit_t>::BasisSymmetricNoNp(int64_t nsites,
                                              Representation const &irrep) try
    : nsites_(nsites), group_action_(irrep.group()), irrep_(irrep),
      raw_ups_size_((int64_t)1 << nsites), raw_dns_size_((int64_t)1 << nsites),
      lintable_ups_(nsites), lintable_dns_(nsites),
      fermi_table_(nsites_, irrep.group()) {
  check_nsites_work_with_bits<bit_t>(nsites_);

  if (nsites < 0) {
    XDIAG_THROW("nsites < 0");
  } else if (nsites != irrep.group().nsites()) {
    XDIAG_THROW("nsites does not match the nsites in PermutationGroup");
  }

  using combinatorics::Subsets;

  std::tie(reps_up_, idces_up_, syms_up_, sym_limits_up_) =
      symmetries::representatives_indices_symmetries_limits<bit_t>(
          combinatorics::SubsetsIndexing<bit_t>(nsites), group_action_);

  if (isreal(irrep)) {
    auto characters = irrep.characters().as<arma::vec>();
    std::tie(dns_storage_, norms_storage_, dns_limits_, ups_offset_, size_) =
        symmetries::electrondns_norms_limits_offset_size(
            reps_up_, Subsets<bit_t>(nsites), group_action_, characters);
  } else {
    auto characters = irrep.characters().as<arma::cx_vec>();
    std::tie(dns_storage_, norms_storage_, dns_limits_, ups_offset_, size_) =
        symmetries::electrondns_norms_limits_offset_size(
            reps_up_, Subsets<bit_t>(nsites), group_action_, characters);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <class bit_t> int64_t BasisSymmetricNoNp<bit_t>::nsites() const {
  return nsites_;
}
template <class bit_t> int64_t BasisSymmetricNoNp<bit_t>::nup() const {
  return nup_;
}
template <class bit_t> int64_t BasisSymmetricNoNp<bit_t>::ndn() const {
  return ndn_;
}

template <class bit_t> int64_t BasisSymmetricNoNp<bit_t>::dim() const {
  return size_;
}
template <class bit_t> int64_t BasisSymmetricNoNp<bit_t>::size() const {
  return size_;
}
template <typename bit_t>
typename BasisSymmetricNoNp<bit_t>::iterator_t
BasisSymmetricNoNp<bit_t>::begin() const {
  return iterator_t(*this, true);
}
template <typename bit_t>
typename BasisSymmetricNoNp<bit_t>::iterator_t
BasisSymmetricNoNp<bit_t>::end() const {
  return iterator_t(*this, false);
}

template <typename bit_t>
int64_t BasisSymmetricNoNp<bit_t>::index(bit_t ups, bit_t dns) const {
  auto syms = syms_ups(ups);
  int64_t idx_ups = index_ups(ups);
  int64_t up_offset = ups_offset(idx_ups);
  // trivial up-stabilizer (likely)
  if (syms.size() == 1) {
    int64_t idx_dns = index_dns(dns);
    return up_offset + idx_dns;
  }
  // non-trivial up-stabilizer (unlikely)
  else {
    auto dnss = dns_for_ups_rep(ups);
    auto [idx_dns, fermi_dn, sym] = index_dns_fermi_sym(dns, syms, dnss);
    return up_offset + idx_dns;
  }
}

template <class bit_t>
GroupActionLookup<bit_t> const &
BasisSymmetricNoNp<bit_t>::group_action() const {
  return group_action_;
}
template <class bit_t>
Representation const &BasisSymmetricNoNp<bit_t>::irrep() const {
  return irrep_;
}

template <typename bit_t>
bool BasisSymmetricNoNp<bit_t>::operator==(
    BasisSymmetricNoNp const &rhs) const {
  return (nsites_ == rhs.nsites_) && (nup_ == rhs.nup_) && (ndn_ == rhs.ndn_) &&
         (group_action_ == rhs.group_action_) && (irrep_ == rhs.irrep_);
}

template <typename bit_t>
bool BasisSymmetricNoNp<bit_t>::operator!=(
    BasisSymmetricNoNp const &rhs) const {
  return !operator==(rhs);
}

template class BasisSymmetricNoNp<uint32_t>;
template class BasisSymmetricNoNp<uint64_t>;

template <typename bit_t>
BasisSymmetricNoNpIterator<bit_t>::BasisSymmetricNoNpIterator(
    BasisSymmetricNoNp<bit_t> const &basis, bool begin)
    : basis_(basis), up_idx_(begin ? 0 : basis.n_rep_ups()), dn_idx_(0) {
  if (basis.dim() == 0) {
    up_idx_ = 0;
  }

  if ((basis.n_rep_ups() > 0) && begin) {
    dns_for_ups_rep_ = basis.dns_for_ups_rep(basis.rep_ups(0));
  }
}

template <typename bit_t>
BasisSymmetricNoNpIterator<bit_t> &
BasisSymmetricNoNpIterator<bit_t>::operator++() {
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
std::pair<bit_t, bit_t> BasisSymmetricNoNpIterator<bit_t>::operator*() const {
  return {basis_.rep_ups(up_idx_), dns_for_ups_rep_[dn_idx_]};
}

template <typename bit_t>
bool BasisSymmetricNoNpIterator<bit_t>::operator!=(
    BasisSymmetricNoNpIterator<bit_t> const &rhs) const {
  return (up_idx_ != rhs.up_idx_) || (dn_idx_ != rhs.dn_idx_);
}

template class BasisSymmetricNoNpIterator<uint32_t>;
template class BasisSymmetricNoNpIterator<uint64_t>;

} // namespace xdiag::basis::electron
