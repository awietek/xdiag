// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "spinhalf_distributed.hpp"

#include <optional>

#include <xdiag/basis/distributed/basis_spinhalf_distributed.hpp>
#include <xdiag/blocks/print_block.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/to_string_generic.hpp>

namespace xdiag {

SpinhalfDistributed::SpinhalfDistributed(int64_t nsites,
                                         RepresentationSet const &irreps) try
    : irreps_(irreps) {
  using namespace basis;
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  }
  if (irreps.group("SitePermutation")) {
    XDIAG_THROW("SpinhalfDistributed does not support permutation symmetries");
  }
  std::optional<int64_t> nup = irreps.charge("nup");
  if (!nup) {
    XDIAG_THROW(
        "SpinhalfDistributed requires a fixed nup (Sz number conservation)");
  }
  if ((*nup < 0) || (*nup > nsites)) {
    XDIAG_THROW("Invalid argument: nup out of range [0, nsites]");
  }

  if (nsites <= 32) {
    basis_ = std::make_shared<BasisSpinhalfDistributed<uint32_t>>(nsites, *nup);
  } else if (nsites <= 64) {
    basis_ = std::make_shared<BasisSpinhalfDistributed<uint64_t>>(nsites, *nup);
  } else {
    XDIAG_THROW("SpinhalfDistributed currently supports up to 64 sites");
  }
}
XDIAG_CATCH

SpinhalfDistributed::SpinhalfDistributed(int64_t nsites, int64_t nup) try
    : SpinhalfDistributed(nsites,
                          RepresentationSet{Representation("nup", nup)}) {}
XDIAG_CATCH

int64_t SpinhalfDistributed::nsites() const { return basis_->nsites(); }
int64_t SpinhalfDistributed::dim() const { return basis_->dim(); }
int64_t SpinhalfDistributed::size() const { return basis_->size(); }
int64_t SpinhalfDistributed::size_max() const { return basis_->size_max(); }
int64_t SpinhalfDistributed::size_min() const { return basis_->size_min(); }
bool SpinhalfDistributed::isreal() const { return irreps_.isreal(); }

bool SpinhalfDistributed::operator==(SpinhalfDistributed const &rhs) const {
  return (irreps_ == rhs.irreps_) && (basis_ == rhs.basis_);
}
bool SpinhalfDistributed::operator!=(SpinhalfDistributed const &rhs) const {
  return !operator==(rhs);
}

SpinhalfDistributedIterator SpinhalfDistributed::begin() const {
  return {this, 0};
}
SpinhalfDistributedIterator SpinhalfDistributed::end() const {
  return {this, size()};
}

RepresentationSet SpinhalfDistributed::irreps() const { return irreps_; }
std::shared_ptr<basis::Basis> const &SpinhalfDistributed::basis() const {
  return basis_;
}

std::ostream &operator<<(std::ostream &out, SpinhalfDistributed const &block) {
  print_block(out, block);
  return out;
}
std::string to_string(SpinhalfDistributed const &block) {
  return to_string_generic(block);
}

SpinhalfDistributedIterator::SpinhalfDistributedIterator(
    SpinhalfDistributed const *block, int64_t idx)
    : idx_(idx) {
  if (idx_ < block->size()) {
    it_ = block->basis()->product_state_iterator();
  }
}

SpinhalfDistributedIterator &SpinhalfDistributedIterator::operator++() {
  it_->advance();
  ++idx_;
  return *this;
}

ProductState SpinhalfDistributedIterator::operator*() const {
  return it_->product_state();
}

bool SpinhalfDistributedIterator::operator==(
    SpinhalfDistributedIterator const &rhs) const {
  return idx_ == rhs.idx_;
}
bool SpinhalfDistributedIterator::operator!=(
    SpinhalfDistributedIterator const &rhs) const {
  return idx_ != rhs.idx_;
}

} // namespace xdiag
