#include "electron.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/utils/logger.h>

namespace hydra {

using namespace indexing;

Electron::Electron(int64_t n_sites)
    : n_sites_(n_sites), charge_conserved_(false), charge_(undefined_qn),
      sz_conserved_(false), sz_(undefined_qn), n_up_(undefined_qn),
      n_dn_(undefined_qn), symmetric_(false), permutation_group_(), irrep_() {

  if (n_sites <= 0) {
    Log.err(
        "Error creating Electron: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ = std::make_shared<ElectronIndexing>(
        electron::IndexingNoNp<uint16_t>(n_sites));
  } else if (n_sites < 32) {
    indexing_ = std::make_shared<ElectronIndexing>(
        electron::IndexingNoNp<uint32_t>(n_sites));
  } else if (n_sites < 64) {
    indexing_ = std::make_shared<ElectronIndexing>(
        electron::IndexingNoNp<uint64_t>(n_sites));
  } else {
    Log.err("Error creating Electron: blocks with more than 64 sites currently "
            "not implemented");
  }
  size_ = indexing::size(*indexing_);
  assert(n_sites >= 0);
}

Electron::Electron(int64_t n_sites, int64_t nup, int64_t ndn)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      symmetric_(false), permutation_group_(), irrep_() {

  if (n_sites <= 0) {
    Log.err(
        "Error creating Electron: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ = std::make_shared<ElectronIndexing>(
        electron::IndexingNp<uint16_t>(n_sites, nup, ndn));
  } else if (n_sites < 32) {
    indexing_ = std::make_shared<ElectronIndexing>(
        electron::IndexingNp<uint32_t>(n_sites, nup, ndn));
  } else if (n_sites < 64) {
    indexing_ = std::make_shared<ElectronIndexing>(
        electron::IndexingNp<uint64_t>(n_sites, nup, ndn));
  } else {
    Log.err("Error creating Electron: blocks with more than 64 sites currently "
            "not implemented");
  }
  size_ = indexing::size(*indexing_);
  utils::check_nup_ndn_electron(n_sites, nup, ndn, "Electron");
}

Electron::Electron(int64_t n_sites, PermutationGroup group, Representation irrep)
    : n_sites_(n_sites), charge_conserved_(false), charge_(undefined_qn),
      sz_conserved_(false), sz_(undefined_qn), n_up_(undefined_qn),
      n_dn_(undefined_qn), symmetric_(true),
      permutation_group_(allowed_subgroup(group, irrep)), irrep_(irrep) {

  if (n_sites <= 0) {
    Log.err(
        "Error creating Electron: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ = std::make_shared<ElectronIndexing>(
        electron::IndexingSymmetricNoNp<uint16_t>(n_sites, group, irrep));
  } else if (n_sites < 32) {
    indexing_ = std::make_shared<ElectronIndexing>(
        electron::IndexingSymmetricNoNp<uint32_t>(n_sites, group, irrep));
  } else if (n_sites < 64) {
    indexing_ = std::make_shared<ElectronIndexing>(
        electron::IndexingSymmetricNoNp<uint64_t>(n_sites, group, irrep));
  } else {
    Log.err("Error creating Electron: blocks with more than 64 sites currently "
            "not implemented");
  }
  size_ = indexing::size(*indexing_);

  utils::check_n_sites(n_sites, group);
}

Electron::Electron(int64_t n_sites, int64_t nup, int64_t ndn, PermutationGroup group,
                   Representation irrep)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      symmetric_(true), permutation_group_(allowed_subgroup(group, irrep)),
      irrep_(irrep) {

  if (n_sites <= 0) {
    Log.err(
        "Error creating Electron: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ = std::make_shared<ElectronIndexing>(
        electron::IndexingSymmetricNp<uint16_t>(n_sites, nup, ndn, group,
                                                irrep));
  } else if (n_sites < 32) {
    indexing_ = std::make_shared<ElectronIndexing>(
        electron::IndexingSymmetricNp<uint32_t>(n_sites, nup, ndn, group,
                                                irrep));
  } else if (n_sites < 64) {
    indexing_ = std::make_shared<ElectronIndexing>(
        electron::IndexingSymmetricNp<uint64_t>(n_sites, nup, ndn, group,
                                                irrep));
  } else {
    Log.err("Error creating Electron: blocks with more than 64 sites currently "
            "not implemented");
  }
  size_ = indexing::size(*indexing_);

  utils::check_nup_ndn_electron(n_sites, nup, ndn, "Electron");
  utils::check_n_sites(n_sites, group);
}

bool Electron::operator==(Electron const &rhs) const {
  return (n_sites_ == rhs.n_sites_) &&
         (charge_conserved_ == rhs.charge_conserved_) &&
         (charge_ == rhs.charge_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}
bool Electron::operator!=(Electron const &rhs) const {
  return !operator==(rhs);
}

indexing::ElectronIndexing const &Electron::indexing() const {
  return *indexing_;
}

} // namespace hydra
