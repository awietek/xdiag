#include "spinhalf.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>
#include <hydra/utils/logger.h>

namespace hydra {

using namespace indexing;

template <typename bit_t>
Spinhalf<bit_t>::Spinhalf(int n_sites)
    : n_sites_(n_sites), sz_conserved_(false), n_up_(undefined_qn),
      n_dn_(undefined_qn), sz_(undefined_qn), symmetric_(false), n_sublat_(0),
      permutation_group_(), irrep_(),
      indexing_(std::make_shared<spinhalf::Indexing<bit_t>>(
          spinhalf::IndexingNoSz<bit_t>(n_sites))),
      size_((idx_t)1 << n_sites) {
  assert(n_sites >= 0);
}

template <typename bit_t>
Spinhalf<bit_t>::Spinhalf(int n_sites, int n_up)
    : n_sites_(n_sites), sz_conserved_(true), n_up_(n_up),
      n_dn_(n_sites - n_up), sz_(n_up_ - n_dn_), symmetric_(false),
      n_sublat_(0), permutation_group_(), irrep_(),
      indexing_(std::make_shared<spinhalf::Indexing<bit_t>>(
          spinhalf::IndexingSz<bit_t>(n_sites, n_up))),
      size_(combinatorics::binomial(n_sites, n_up)) {
  assert(n_sites >= 0);
  utils::check_nup_spinhalf(n_sites, n_up, "Spinhalf");
}

template <typename bit_t>
Spinhalf<bit_t>::Spinhalf(int n_sites, PermutationGroup permutation_group,
                          Representation irrep, int n_sublat)
    : n_sites_(n_sites), sz_conserved_(false), n_up_(undefined_qn),
      n_dn_(undefined_qn), sz_(undefined_qn), symmetric_(true),
      n_sublat_(n_sublat),
      permutation_group_(allowed_subgroup(permutation_group, irrep)),
      irrep_(irrep) {

  utils::check_n_sites(n_sites, permutation_group);
  if (n_sublat == 0) {
    indexing_ = std::make_shared<spinhalf::Indexing<bit_t>>(
        spinhalf::IndexingSymmetricNoSz<bit_t>(n_sites, permutation_group,
                                               irrep));
  } else if (n_sublat == 1) {
    indexing_ = std::make_shared<spinhalf::Indexing<bit_t>>(
        spinhalf::IndexingSublattice<bit_t, 1>(n_sites, permutation_group,
                                               irrep));
  } else if (n_sublat == 2) {
    indexing_ = std::make_shared<spinhalf::Indexing<bit_t>>(
        spinhalf::IndexingSublattice<bit_t, 2>(n_sites, permutation_group,
                                               irrep));
  } else if (n_sublat == 3) {
    indexing_ = std::make_shared<spinhalf::Indexing<bit_t>>(
        spinhalf::IndexingSublattice<bit_t, 3>(n_sites, permutation_group,
                                               irrep));
  } else if (n_sublat == 4) {
    indexing_ = std::make_shared<spinhalf::Indexing<bit_t>>(
        spinhalf::IndexingSublattice<bit_t, 4>(n_sites, permutation_group,
                                               irrep));
  } else if (n_sublat == 5) {
    indexing_ = std::make_shared<spinhalf::Indexing<bit_t>>(
        spinhalf::IndexingSublattice<bit_t, 5>(n_sites, permutation_group,
                                               irrep));
  } else {
    Log.err("Error creating Spinhalf: invalid n_sublat specified. Must be "
            "eiter 0 (so sublattice coding) or between 1 and 5. Got {}",
            n_sublat);
  }
  size_ = spinhalf::size(*indexing_);
}

template <typename bit_t>
Spinhalf<bit_t>::Spinhalf(int n_sites, int n_up,
                          PermutationGroup permutation_group,
                          Representation irrep, int n_sublat)
    : n_sites_(n_sites), sz_conserved_(true), n_up_(n_up),
      n_dn_(n_sites - n_up), sz_(n_up_ - n_dn_), symmetric_(true),
      n_sublat_(n_sublat),
      permutation_group_(allowed_subgroup(permutation_group, irrep)),
      irrep_(irrep) {
  utils::check_nup_spinhalf(n_sites, n_up, "Spinhalf");
  utils::check_n_sites(n_sites, permutation_group);

  if (n_sublat == 0) {
    indexing_ = std::make_shared<spinhalf::Indexing<bit_t>>(
        spinhalf::IndexingSymmetricSz<bit_t>(n_sites, n_up, permutation_group,
                                             irrep));
  } else if (n_sublat == 1) {
    indexing_ = std::make_shared<spinhalf::Indexing<bit_t>>(
        spinhalf::IndexingSublattice<bit_t, 1>(n_sites, n_up, permutation_group,
                                               irrep));
  } else if (n_sublat == 2) {
    indexing_ = std::make_shared<spinhalf::Indexing<bit_t>>(
        spinhalf::IndexingSublattice<bit_t, 2>(n_sites, n_up, permutation_group,
                                               irrep));
  } else if (n_sublat == 3) {
    indexing_ = std::make_shared<spinhalf::Indexing<bit_t>>(
        spinhalf::IndexingSublattice<bit_t, 3>(n_sites, n_up, permutation_group,
                                               irrep));
  } else if (n_sublat == 4) {
    indexing_ = std::make_shared<spinhalf::Indexing<bit_t>>(
        spinhalf::IndexingSublattice<bit_t, 4>(n_sites, n_up, permutation_group,
                                               irrep));
  } else if (n_sublat == 5) {
    indexing_ = std::make_shared<spinhalf::Indexing<bit_t>>(
        spinhalf::IndexingSublattice<bit_t, 5>(n_sites, n_up, permutation_group,
                                               irrep));
  } else {
    Log.err("Error creating Spinhalf: invalid n_sublat specified. Must be "
            "eiter 0 (so sublattice coding) or between 1 and 5. Got {}",
            n_sublat);
  }
  size_ = spinhalf::size(*indexing_);
}

template <typename bit_t>
bool Spinhalf<bit_t>::operator==(Spinhalf<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (symmetric_ == rhs.symmetric_) && (n_sublat_ == rhs.n_sublat_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}

template <typename bit_t>
bool Spinhalf<bit_t>::operator!=(Spinhalf<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template <typename bit_t>
indexing::spinhalf::Indexing<bit_t> const &Spinhalf<bit_t>::indexing() const {
  return *indexing_;
}

template class Spinhalf<uint16_t>;
template class Spinhalf<uint32_t>;
template class Spinhalf<uint64_t>;

} // namespace hydra
