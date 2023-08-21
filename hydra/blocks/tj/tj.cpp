#include "tj.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/utils/logger.h>

namespace hydra {

using namespace indexing;

tJ::tJ(int64_t n_sites, int64_t nup, int64_t ndn)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      symmetric_(false), permutation_group_(), irrep_() {

  if (n_sites <= 0) {
    Log.err("Error creating tJ: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ = std::make_shared<tJIndexing>(
        tj::IndexingNp<uint16_t>(n_sites, nup, ndn));
  } else if (n_sites < 32) {
    indexing_ = std::make_shared<tJIndexing>(
        tj::IndexingNp<uint32_t>(n_sites, nup, ndn));
  } else if (n_sites < 64) {
    indexing_ = std::make_shared<tJIndexing>(
        tj::IndexingNp<uint64_t>(n_sites, nup, ndn));
  } else {
    Log.err("Error creating tJ: blocks with more than 64 sites currently "
            "not implemented");
  }
  size_ = indexing::size(*indexing_);
  utils::check_nup_ndn_tj(n_sites, nup, ndn, "tJ");
}

tJ::tJ(int64_t n_sites, int64_t nup, int64_t ndn, PermutationGroup group,
       Representation irrep)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      symmetric_(true), permutation_group_(allowed_subgroup(group, irrep)),
      irrep_(irrep) {

  if (n_sites <= 0) {
    Log.err("Error creating tJ: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ = std::make_shared<tJIndexing>(
        tj::IndexingSymmetricNp<uint16_t>(n_sites, nup, ndn, group, irrep));
  } else if (n_sites < 32) {
    indexing_ = std::make_shared<tJIndexing>(
        tj::IndexingSymmetricNp<uint32_t>(n_sites, nup, ndn, group, irrep));
  } else if (n_sites < 64) {
    indexing_ = std::make_shared<tJIndexing>(
        tj::IndexingSymmetricNp<uint64_t>(n_sites, nup, ndn, group, irrep));
  } else {
    Log.err("Error creating tJ: blocks with more than 64 sites currently "
            "not implemented");
  }
  size_ = indexing::size(*indexing_);
  utils::check_nup_ndn_tj(n_sites, nup, ndn, "tJ");
  utils::check_n_sites(n_sites, group);
}

bool tJ::operator==(tJ const &rhs) const {
  return (n_sites_ == rhs.n_sites_) &&
         (charge_conserved_ == rhs.charge_conserved_) &&
         (charge_ == rhs.charge_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}
bool tJ::operator!=(tJ const &rhs) const { return !operator==(rhs); }

indexing::tJIndexing const &tJ::indexing() const { return *indexing_; }

} // namespace hydra
