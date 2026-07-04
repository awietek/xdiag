// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "electron_distributed.hpp"

#include <optional>

#include <xdiag/basis/distributed/basis_electron_distributed.hpp>
#include <xdiag/blocks/print_block.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/to_string_generic.hpp>

namespace xdiag {

ElectronDistributed::ElectronDistributed(int64_t nsites,
                                         RepresentationSet const &irreps) try
    : irreps_(irreps) {
  using namespace basis;
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  }
  if (irreps.group("SitePermutation")) {
    XDIAG_THROW("ElectronDistributed does not support permutation symmetries");
  }
  std::optional<int64_t> nup = irreps.charge("nup");
  std::optional<int64_t> ndn = irreps.charge("ndn");
  if (!nup || !ndn) {
    XDIAG_THROW("ElectronDistributed requires fixed nup and ndn (number "
                "conservation in both spin species)");
  }
  if ((*nup < 0) || (*nup > nsites) || (*ndn < 0) || (*ndn > nsites)) {
    XDIAG_THROW("Invalid argument: nup or ndn out of range [0, nsites]");
  }

  if (nsites <= 32) {
    basis_ =
        std::make_shared<BasisElectronDistributed<uint32_t>>(nsites, *nup, *ndn);
  } else if (nsites <= 64) {
    basis_ =
        std::make_shared<BasisElectronDistributed<uint64_t>>(nsites, *nup, *ndn);
  } else {
    XDIAG_THROW("ElectronDistributed currently supports up to 64 sites");
  }
}
XDIAG_CATCH

ElectronDistributed::ElectronDistributed(int64_t nsites, int64_t nup,
                                         int64_t ndn) try
    : ElectronDistributed(nsites,
                          RepresentationSet{Representation("nup", nup),
                                            Representation("ndn", ndn)}) {}
XDIAG_CATCH

int64_t ElectronDistributed::nsites() const { return basis_->nsites(); }
int64_t ElectronDistributed::dim() const { return basis_->dim(); }
int64_t ElectronDistributed::size() const { return basis_->size(); }
int64_t ElectronDistributed::size_max() const { return basis_->size_max(); }
int64_t ElectronDistributed::size_min() const { return basis_->size_min(); }
bool ElectronDistributed::isreal() const { return irreps_.isreal(); }
int64_t ElectronDistributed::index(ProductState const &pstate) const {
  return basis_->index(pstate);
}
  
bool ElectronDistributed::operator==(ElectronDistributed const &rhs) const {
  return (irreps_ == rhs.irreps_) && (basis_ == rhs.basis_);
}
bool ElectronDistributed::operator!=(ElectronDistributed const &rhs) const {
  return !operator==(rhs);
}

ElectronDistributedIterator ElectronDistributed::begin() const {
  return {this, 0};
}
ElectronDistributedIterator ElectronDistributed::end() const {
  return {this, size()};
}

RepresentationSet ElectronDistributed::irreps() const { return irreps_; }
std::shared_ptr<basis::Basis> const &ElectronDistributed::basis() const {
  return basis_;
}

std::ostream &operator<<(std::ostream &out, ElectronDistributed const &block) {
  print_block(out, block);
  return out;
}
std::string to_string(ElectronDistributed const &block) {
  return to_string_generic(block);
}

ElectronDistributedIterator::ElectronDistributedIterator(
    ElectronDistributed const *block, int64_t idx)
    : idx_(idx) {
  if (idx_ < block->size()) {
    it_ = block->basis()->product_state_iterator();
  }
}

ElectronDistributedIterator &ElectronDistributedIterator::operator++() {
  it_->advance();
  ++idx_;
  return *this;
}

ProductState ElectronDistributedIterator::operator*() const {
  return it_->product_state();
}

bool ElectronDistributedIterator::operator==(
    ElectronDistributedIterator const &rhs) const {
  return idx_ == rhs.idx_;
}
bool ElectronDistributedIterator::operator!=(
    ElectronDistributedIterator const &rhs) const {
  return idx_ != rhs.idx_;
}

} // namespace xdiag
