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

} // namespace xdiag::basis::electron
