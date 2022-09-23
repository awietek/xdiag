#include "tj.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/utils/logger.h>

namespace hydra {

using namespace indexing;

template <class bit_t>
tJ<bit_t>::tJ(int n_sites, int nup, int ndn)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      symmetric_(false), permutation_group_(), irrep_(),
      indexing_(std::make_shared<tj::Indexing<bit_t>>(
          tj::IndexingNp<bit_t>(n_sites, nup, ndn))),
      size_(tj::size(*indexing_)) {
  assert(n_sites >= 0);
  utils::check_nup_ndn_tj(n_sites, nup, ndn, "tJ");
}

template <typename bit_t>
tJ<bit_t>::tJ(int n_sites, int nup, int ndn, PermutationGroup permutation_group,
              Representation irrep)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      symmetric_(true),
      permutation_group_(allowed_subgroup(permutation_group, irrep)),
      irrep_(irrep), indexing_(std::make_shared<tj::Indexing<bit_t>>(
                         tj::IndexingNp<bit_t>(n_sites, nup, ndn))),
      size_(tj::size(*indexing_)) {
  utils::check_nup_ndn_tj(n_sites, nup, ndn, "tJ");
  utils::check_n_sites(n_sites, permutation_group);
}

template <typename bit_t>
bool tJ<bit_t>::operator==(tJ<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) &&
         (charge_conserved_ == rhs.charge_conserved_) &&
         (charge_ == rhs.charge_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}
template <class bit_t> bool tJ<bit_t>::operator!=(tJ<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template <typename bit_t>
indexing::tj::Indexing<bit_t> const &tJ<bit_t>::indexing() const {
  return *indexing_;
}

template class tJ<uint16_t>;
template class tJ<uint32_t>;
template class tJ<uint64_t>;

} // namespace hydra
