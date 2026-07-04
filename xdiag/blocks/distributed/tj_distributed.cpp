// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "tj_distributed.hpp"

#include <optional>

#include <xdiag/basis/distributed/basis_tj_distributed.hpp>
#include <xdiag/blocks/print_block.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/to_string_generic.hpp>

namespace xdiag {

tJDistributed::tJDistributed(int64_t nsites,
                             RepresentationSet const &irreps) try
    : irreps_(irreps) {
  using namespace basis;
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  }
  if (irreps.group("SitePermutation")) {
    XDIAG_THROW("tJDistributed does not support permutation symmetries");
  }
  std::optional<int64_t> nup = irreps.charge("nup");
  std::optional<int64_t> ndn = irreps.charge("ndn");
  if (!nup || !ndn) {
    XDIAG_THROW("tJDistributed requires fixed nup and ndn (number conservation "
                "in both spin species)");
  }
  if ((*nup < 0) || (*ndn < 0)) {
    XDIAG_THROW("Invalid argument: nup < 0 or ndn < 0");
  }
  if ((*nup + *ndn) > nsites) {
    XDIAG_THROW("Invalid argument: nup + ndn > nsites (no double occupancy)");
  }

  if (nsites <= 32) {
    basis_ = std::make_shared<BasistJDistributed<uint32_t>>(nsites, *nup, *ndn);
  } else if (nsites <= 64) {
    basis_ = std::make_shared<BasistJDistributed<uint64_t>>(nsites, *nup, *ndn);
  } else {
    XDIAG_THROW("tJDistributed currently supports up to 64 sites");
  }
}
XDIAG_CATCH

tJDistributed::tJDistributed(int64_t nsites, int64_t nup, int64_t ndn) try
    : tJDistributed(nsites, RepresentationSet{Representation("nup", nup),
                                              Representation("ndn", ndn)}) {}
XDIAG_CATCH

int64_t tJDistributed::nsites() const { return basis_->nsites(); }
int64_t tJDistributed::dim() const { return basis_->dim(); }
int64_t tJDistributed::size() const { return basis_->size(); }
int64_t tJDistributed::size_max() const { return basis_->size_max(); }
int64_t tJDistributed::size_min() const { return basis_->size_min(); }
bool tJDistributed::isreal() const { return irreps_.isreal(); }
int64_t tJDistributed::index(ProductState const &pstate) const {
  return basis_->index(pstate);
}

bool tJDistributed::operator==(tJDistributed const &rhs) const {
  return (irreps_ == rhs.irreps_) && (basis_ == rhs.basis_);
}
bool tJDistributed::operator!=(tJDistributed const &rhs) const {
  return !operator==(rhs);
}

tJDistributedIterator tJDistributed::begin() const { return {this, 0}; }
tJDistributedIterator tJDistributed::end() const { return {this, size()}; }

RepresentationSet tJDistributed::irreps() const { return irreps_; }
std::shared_ptr<basis::Basis> const &tJDistributed::basis() const {
  return basis_;
}

std::ostream &operator<<(std::ostream &out, tJDistributed const &block) {
  print_block(out, block);
  return out;
}
std::string to_string(tJDistributed const &block) {
  return to_string_generic(block);
}

tJDistributedIterator::tJDistributedIterator(tJDistributed const *block,
                                             int64_t idx)
    : idx_(idx) {
  if (idx_ < block->size()) {
    it_ = block->basis()->product_state_iterator();
  }
}

tJDistributedIterator &tJDistributedIterator::operator++() {
  it_->advance();
  ++idx_;
  return *this;
}

ProductState tJDistributedIterator::operator*() const {
  return it_->product_state();
}

bool tJDistributedIterator::operator==(tJDistributedIterator const &rhs) const {
  return idx_ == rhs.idx_;
}
bool tJDistributedIterator::operator!=(tJDistributedIterator const &rhs) const {
  return idx_ != rhs.idx_;
}

} // namespace xdiag
