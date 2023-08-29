#include "basis_symmetric_no_np.h"

#include <hydra/combinatorics/subsets.h>
#include <hydra/symmetries/operations/group_action_operations.h>
#include <hydra/symmetries/operations/representative_list.h>
#include <hydra/symmetries/operations/symmetry_operations.h>

namespace hydra::basis::electron {

template <class bit_t>
BasisSymmetricNoNp<bit_t>::BasisSymmetricNoNp(int64_t n_sites,
                                              PermutationGroup group,
                                              Representation irrep)
    : n_sites_(n_sites), group_action_(allowed_subgroup(group, irrep)),
      irrep_(irrep), raw_ups_size_((idx_t)1 << n_sites),
      raw_dns_size_((idx_t)1 << n_sites), lintable_ups_(n_sites),
      lintable_dns_(n_sites), fermi_table_(n_sites_, group) {
  if (n_sites < 0) {
    throw(std::invalid_argument("n_sites < 0"));
  } else if (n_sites != group.n_sites()) {
    throw(std::logic_error(
        "n_sites does not match the n_sites in PermutationGroup"));
  }

  using combinatorics::Subsets;

  std::tie(reps_up_, idces_up_, syms_up_, sym_limits_up_) =
      symmetries::representatives_indices_symmetries_limits<bit_t>(
          combinatorics::SubsetsIndexing<bit_t>(n_sites), group_action_);

  std::tie(dns_storage_, norms_storage_, dns_limits_, ups_offset_, size_) =
      symmetries::electron_dns_norms_limits_offset_size(
          reps_up_, Subsets<bit_t>(n_sites), group_action_, irrep_);
}

template class BasisSymmetricNoNp<uint16_t>;
template class BasisSymmetricNoNp<uint32_t>;
template class BasisSymmetricNoNp<uint64_t>;

} // namespace hydra::basis::electron
