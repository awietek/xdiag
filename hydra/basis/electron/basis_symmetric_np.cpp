#include "basis_symmetric_np.h"

#include <hydra/combinatorics/combinations_indexing.h>
#include <hydra/combinatorics/fermi_table.h>

#include <hydra/symmetries/operations/group_action_operations.h>
#include <hydra/symmetries/operations/representative_list.h>
#include <hydra/symmetries/operations/symmetry_operations.h>

namespace hydra::basis::electron {

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

template class BasisSymmetricNp<uint16_t>;
template class BasisSymmetricNp<uint32_t>;
template class BasisSymmetricNp<uint64_t>;

} // namespace hydra::basis::electron
