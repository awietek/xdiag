#include "spinhalf.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>
#include <hydra/utils/logger.h>

namespace hydra {

using namespace indexing;

Spinhalf::Spinhalf(int64_t n_sites)
    : n_sites_(n_sites), sz_conserved_(false), n_up_(undefined_qn),
      n_dn_(undefined_qn), sz_(undefined_qn), symmetric_(false), n_sublat_(0),
      permutation_group_(), irrep_(), size_((idx_t)1 << n_sites) {

  if (n_sites <= 0) {
    Log.err(
        "Error creating Spinhalf: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingNoSz<uint16_t>(n_sites));
  } else if (n_sites < 32) {
    indexing_ = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingNoSz<uint32_t>(n_sites));
  } else if (n_sites < 64) {
    indexing_ = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingNoSz<uint64_t>(n_sites));
  } else {
    Log.err("Error creating Spinhalf: blocks with more than 64 sites currently "
            "not implemented");
  }
}

Spinhalf::Spinhalf(int64_t n_sites, int64_t n_up)
    : n_sites_(n_sites), sz_conserved_(true), n_up_(n_up),
      n_dn_(n_sites - n_up), sz_(n_up_ - n_dn_), symmetric_(false),
      n_sublat_(0), permutation_group_(), irrep_(),
      size_(combinatorics::binomial(n_sites, n_up)) {

  if (n_sites <= 0) {
    Log.err(
        "Error creating Spinhalf: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSz<uint16_t>(n_sites, n_up));
  } else if (n_sites < 32) {
    indexing_ = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSz<uint32_t>(n_sites, n_up));
  } else if (n_sites < 64) {
    indexing_ = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSz<uint64_t>(n_sites, n_up));
  } else {
    Log.err("Error creating Spinhalf: blocks with more than 64 sites currently "
            "not implemented");
  }
  utils::check_nup_spinhalf(n_sites, n_up, "Spinhalf");
}

template <typename bit_t>
std::shared_ptr<SpinhalfIndexing>
make_spinhalf_indexing_no_sz(int64_t n_sites, PermutationGroup const &group,
                             Representation const &irrep, int64_t n_sublat) {
  std::shared_ptr<SpinhalfIndexing> indexing;
  if (n_sublat == 0) {
    indexing = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSymmetricNoSz<bit_t>(n_sites, group, irrep));
  } else if (n_sublat == 1) {
    indexing = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSublattice<bit_t, 1>(n_sites, group, irrep));
  } else if (n_sublat == 2) {
    indexing = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSublattice<bit_t, 2>(n_sites, group, irrep));
  } else if (n_sublat == 3) {
    indexing = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSublattice<bit_t, 3>(n_sites, group, irrep));
  } else if (n_sublat == 4) {
    indexing = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSublattice<bit_t, 4>(n_sites, group, irrep));
  } else if (n_sublat == 5) {
    indexing = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSublattice<bit_t, 5>(n_sites, group, irrep));
  } else {
    Log.err("Error creating Spinhalf: invalid n_sublat specified. Must be "
            "eiter 0 (so sublattice coding) or between 1 and 5. Got {}",
            n_sublat);
  }
  return indexing;
}

Spinhalf::Spinhalf(int64_t n_sites, PermutationGroup group, Representation irrep,
                   int64_t n_sublat)
    : n_sites_(n_sites), sz_conserved_(false), n_up_(undefined_qn),
      n_dn_(undefined_qn), sz_(undefined_qn), symmetric_(true),
      n_sublat_(n_sublat), permutation_group_(allowed_subgroup(group, irrep)),
      irrep_(irrep) {

  utils::check_n_sites(n_sites, group);

  if (n_sites <= 0) {
    Log.err(
        "Error creating Spinhalf: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ =
        make_spinhalf_indexing_no_sz<uint16_t>(n_sites, group, irrep, n_sublat);
  } else if (n_sites < 32) {
    indexing_ =
        make_spinhalf_indexing_no_sz<uint32_t>(n_sites, group, irrep, n_sublat);
  } else if (n_sites < 64) {
    indexing_ =
        make_spinhalf_indexing_no_sz<uint64_t>(n_sites, group, irrep, n_sublat);
  } else {
    Log.err("Error creating Spinhalf: blocks with more than 64 sites currently "
            "not implemented");
  }
  size_ = indexing::size(*indexing_);
}

template <typename bit_t>
std::shared_ptr<SpinhalfIndexing>
make_spinhalf_indexing_sz(int64_t n_sites, int64_t n_up, PermutationGroup const &group,
                          Representation const &irrep, int64_t n_sublat) {
  std::shared_ptr<SpinhalfIndexing> indexing;
  if (n_sublat == 0) {
    indexing = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSymmetricSz<bit_t>(n_sites, n_up, group, irrep));
  } else if (n_sublat == 1) {
    indexing = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSublattice<bit_t, 1>(n_sites, n_up, group, irrep));
  } else if (n_sublat == 2) {
    indexing = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSublattice<bit_t, 2>(n_sites, n_up, group, irrep));
  } else if (n_sublat == 3) {
    indexing = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSublattice<bit_t, 3>(n_sites, n_up, group, irrep));
  } else if (n_sublat == 4) {
    indexing = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSublattice<bit_t, 4>(n_sites, n_up, group, irrep));
  } else if (n_sublat == 5) {
    indexing = std::make_shared<SpinhalfIndexing>(
        spinhalf::IndexingSublattice<bit_t, 5>(n_sites, n_up, group, irrep));
  } else {
    Log.err("Error creating Spinhalf: invalid n_sublat specified. Must be "
            "eiter 0 (so sublattice coding) or between 1 and 5. Got {}",
            n_sublat);
  }
  return indexing;
}

Spinhalf::Spinhalf(int64_t n_sites, int64_t n_up, PermutationGroup group,
                   Representation irrep, int64_t n_sublat)
    : n_sites_(n_sites), sz_conserved_(true), n_up_(n_up),
      n_dn_(n_sites - n_up), sz_(n_up_ - n_dn_), symmetric_(true),
      n_sublat_(n_sublat), permutation_group_(allowed_subgroup(group, irrep)),
      irrep_(irrep) {
  utils::check_nup_spinhalf(n_sites, n_up, "Spinhalf");
  utils::check_n_sites(n_sites, group);

  if (n_sites <= 0) {
    Log.err(
        "Error creating Spinhalf: number of sites must be a positive integer");
  } else if (n_sites < 16) {
    indexing_ = make_spinhalf_indexing_sz<uint16_t>(n_sites, n_up, group, irrep,
                                                    n_sublat);
  } else if (n_sites < 32) {
    indexing_ = make_spinhalf_indexing_sz<uint32_t>(n_sites, n_up, group, irrep,
                                                    n_sublat);
  } else if (n_sites < 64) {
    indexing_ = make_spinhalf_indexing_sz<uint64_t>(n_sites, n_up, group, irrep,
                                                    n_sublat);
  } else {
    Log.err("Error creating Spinhalf: blocks with more than 64 sites currently "
            "not implemented");
  }

  size_ = indexing::size(*indexing_);
}

bool Spinhalf::operator==(Spinhalf const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (symmetric_ == rhs.symmetric_) && (n_sublat_ == rhs.n_sublat_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}

bool Spinhalf::operator!=(Spinhalf const &rhs) const {
  return !operator==(rhs);
}

indexing::SpinhalfIndexing const &Spinhalf::indexing() const {
  return *indexing_;
}

} // namespace hydra
