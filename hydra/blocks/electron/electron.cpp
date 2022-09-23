#include "electron.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/utils/logger.h>

namespace hydra {

using namespace indexing;

template <class bit_t>
Electron<bit_t>::Electron(int n_sites)
    : n_sites_(n_sites), charge_conserved_(false), charge_(undefined_qn),
      sz_conserved_(false), sz_(undefined_qn), n_up_(undefined_qn),
      n_dn_(undefined_qn), symmetric_(false), permutation_group_(), irrep_(),
      indexing_(std::make_shared<electron::Indexing<bit_t>>(
          electron::IndexingNoNp<bit_t>(n_sites))),
      size_(electron::size(*indexing_)) {
  assert(n_sites >= 0);
}

template <class bit_t>
Electron<bit_t>::Electron(int n_sites, int nup, int ndn)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      symmetric_(false), permutation_group_(), irrep_(),
      indexing_(std::make_shared<electron::Indexing<bit_t>>(
          electron::IndexingNp<bit_t>(n_sites, nup, ndn))),
      size_(electron::size(*indexing_)) {
  assert(n_sites >= 0);
  utils::check_nup_ndn_electron(n_sites, nup, ndn, "Electron");
}

template <typename bit_t>
Electron<bit_t>::Electron(int n_sites, PermutationGroup permutation_group,
                          Representation irrep)
    : n_sites_(n_sites), charge_conserved_(false), charge_(undefined_qn),
      sz_conserved_(false), sz_(undefined_qn), n_up_(undefined_qn),
      n_dn_(undefined_qn), symmetric_(true),
      permutation_group_(allowed_subgroup(permutation_group, irrep)),
      irrep_(irrep), indexing_(std::make_shared<electron::Indexing<bit_t>>(
                         electron::IndexingSymmetricNoNp<bit_t>(
                             n_sites, permutation_group, irrep))),
      size_(electron::size(*indexing_)) {
  utils::check_n_sites(n_sites, permutation_group);
}

template <typename bit_t>
Electron<bit_t>::Electron(int n_sites, int nup, int ndn,
                          PermutationGroup permutation_group,
                          Representation irrep)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      symmetric_(true),
      permutation_group_(allowed_subgroup(permutation_group, irrep)),
      irrep_(irrep), indexing_(std::make_shared<electron::Indexing<bit_t>>(
                         electron::IndexingSymmetricNp<bit_t>(
                             n_sites, nup, ndn, permutation_group, irrep))),
      size_(electron::size(*indexing_)) {
  utils::check_nup_ndn_electron(n_sites, nup, ndn, "Electron");
  utils::check_n_sites(n_sites, permutation_group);
}

template <typename bit_t>
bool Electron<bit_t>::operator==(Electron<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) &&
         (charge_conserved_ == rhs.charge_conserved_) &&
         (charge_ == rhs.charge_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}
template <class bit_t>
bool Electron<bit_t>::operator!=(Electron<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template <typename bit_t>
indexing::electron::Indexing<bit_t> const &Electron<bit_t>::indexing() const {
  return *indexing_;
}

template class Electron<uint16_t>;
template class Electron<uint32_t>;
template class Electron<uint64_t>;

} // namespace hydra
