#include "basis_symmetric_np.hpp"

#include <xdiag/combinatorics/combinations_indexing.hpp>
#include <xdiag/combinatorics/fermi_table.hpp>

#include <xdiag/symmetries/operations/group_action_operations.hpp>
#include <xdiag/symmetries/operations/representative_list.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>

namespace xdiag::basis::electron {

template <class bit_t>
BasisSymmetricNp<bit_t>::BasisSymmetricNp(int64_t n_sites, int64_t nup,
                                          int64_t ndn,
                                          Representation const &irrep) try
    : n_sites_(n_sites), n_up_(nup), n_dn_(ndn), group_action_(irrep.group()),
      irrep_(irrep), raw_ups_size_(combinatorics::binomial(n_sites, nup)),
      raw_dns_size_(combinatorics::binomial(n_sites, ndn)),
      lintable_ups_(n_sites, nup), lintable_dns_(n_sites, ndn),
      fermi_table_ups_(n_sites, nup, irrep.group()),
      fermi_table_dns_(n_sites, ndn, irrep.group()) {
  check_n_sites_work_with_bits<bit_t>(n_sites_);

  using combinatorics::Combinations;
  if (n_sites < 0) {
    XDIAG_THROW("n_sites < 0");
  } else if ((nup < 0) || (ndn < 0)) {
    XDIAG_THROW("nup < 0 or ndn < 0");
  } else if (n_sites < 0) {
    XDIAG_THROW("n_sites < 0");
  } else if (n_sites != irrep.group().n_sites()) {
    XDIAG_THROW("n_sites does not match the n_sites in PermutationGroup");
  }

  std::tie(reps_up_, idces_up_, syms_up_, sym_limits_up_) =
      symmetries::representatives_indices_symmetries_limits<bit_t>(
          combinatorics::CombinationsIndexing<bit_t>(n_sites, nup),
          group_action_);

  if (isreal(irrep)) {
    auto characters = irrep.characters().as<arma::vec>();
    std::tie(dns_storage_, norms_storage_, dns_limits_, ups_offset_, size_) =
        symmetries::electron_dns_norms_limits_offset_size(
            reps_up_, Combinations<bit_t>(n_sites, ndn), group_action_,
            characters);
  } else {
    auto characters = irrep.characters().as<arma::cx_vec>();
    std::tie(dns_storage_, norms_storage_, dns_limits_, ups_offset_, size_) =
        symmetries::electron_dns_norms_limits_offset_size(
            reps_up_, Combinations<bit_t>(n_sites, ndn), group_action_,
            characters);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
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
template <typename bit_t>
int64_t BasisSymmetricNp<bit_t>::index(bit_t ups, bit_t dns) const {
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
