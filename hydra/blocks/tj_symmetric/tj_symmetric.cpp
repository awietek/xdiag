#include "tj_symmetric.h"

#include <cassert>
#include <tuple>

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/symmetries/fermi_sign.h>
#include <hydra/symmetries/symmetry_operations.h>

namespace hydra {

template <typename bit_t>
tJSymmetric<bit_t>::tJSymmetric(int n_sites, int nup, int ndn,
                                PermutationGroup permutation_group,
                                Representation irrep)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      permutation_group_(
          (irrep.allowed_symmetries().size() > 0)
              ? permutation_group.subgroup(irrep.allowed_symmetries())
              : permutation_group),
      irrep_(irrep), indexing_(n_sites, nup, ndn, permutation_group_, irrep_) {
  utils::check_nup_ndn_tj(n_sites, nup, ndn, "tJSymmetric");
  size_ = indexing_.size();
}

template <typename bit_t>
bool tJSymmetric<bit_t>::operator==(tJSymmetric<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) &&
         (charge_conserved_ == rhs.charge_conserved_) &&
         (charge_ == rhs.charge_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}

template <typename bit_t>
bool tJSymmetric<bit_t>::operator!=(tJSymmetric<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class tJSymmetric<uint16_t>;
template class tJSymmetric<uint32_t>;
template class tJSymmetric<uint64_t>;

} // namespace hydra
