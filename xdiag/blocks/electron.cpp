// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "electron.hpp"

#include <regex>

#include <xdiag/basis/basis_electron.hpp>
#include <xdiag/basis/basis_electron_symmetric.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/blocks/print_block.hpp>
#include <xdiag/random/hash.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>

namespace xdiag {

// Builds the up and dn enumerations (same type) and calls f(enum_up, enum_dn).
// Number conservation requires both nup and ndn to be fixed, or neither.
template <typename F>
static void dispatch_enumeration(int64_t nsites, std::optional<int64_t> nup,
                                 std::optional<int64_t> ndn, F &&f) {
  using namespace bits;
  using namespace combinatorics;
  if (nup && ndn) { // fixed nup and ndn; LinTable for fast lookup up to 42 sites
    int64_t n = nsites;
    int64_t ku = *nup;
    int64_t kd = *ndn;
    if (nsites <= 32) {
      f(LinTable<uint32_t>(n, ku), LinTable<uint32_t>(n, kd));
    } else if (nsites <= 42) {
      f(LinTable<uint64_t>(n, ku), LinTable<uint64_t>(n, kd));
    } else if (nsites <= 64) {
      f(Combinations<uint64_t>(n, ku), Combinations<uint64_t>(n, kd));
    } else {
      f(Combinations<BitsetDynamic>(n, ku), Combinations<BitsetDynamic>(n, kd));
    }
  } else if (!nup && !ndn) { // no number conservation
    if (nsites <= 32) {
      f(Subsets<uint32_t>(nsites), Subsets<uint32_t>(nsites));
    } else if (nsites <= 64) {
      f(Subsets<uint64_t>(nsites), Subsets<uint64_t>(nsites));
    } else {
      XDIAG_THROW(
          "Unsupported nsites > 64 for non-number conserving Electron block");
    }
  } else {
    XDIAG_THROW(
        "Electron block requires both nup and ndn to be set, or neither");
  }
}

Electron::Electron(int64_t nsites, RepresentationSet const &irreps) try
    : irreps_(irreps) {
  using namespace basis;

  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  }

  std::optional<int64_t> nup = irreps.charge("nup");
  std::optional<int64_t> ndn = irreps.charge("ndn");
  if (nup && ((*nup < 0) || (*nup > nsites))) {
    XDIAG_THROW("Invalid argument: nup out of range [0, nsites]");
  }
  if (ndn && ((*ndn < 0) || (*ndn > nsites))) {
    XDIAG_THROW("Invalid argument: ndn out of range [0, nsites]");
  }

  std::optional<PermutationGroup> group = irreps.group("SitePermutation");
  std::optional<Vector> characters = irreps.characters("SitePermutation");

  if (group) { // permutation symmetry -> BasisElectronSymmetric
    dispatch_enumeration(nsites, nup, ndn, [&](auto &&enum_up, auto &&enum_dn) {
      using enum_t = std::decay_t<decltype(enum_up)>;
      basis_ = std::make_shared<BasisElectronSymmetric<enum_t>>(
          enum_up, enum_dn, *group, *characters);
    });
  } else { // no permutation symmetry -> BasisElectron
    dispatch_enumeration(nsites, nup, ndn, [&](auto &&enum_up, auto &&enum_dn) {
      using enum_t = std::decay_t<decltype(enum_up)>;
      basis_ = std::make_shared<BasisElectron<enum_t>>(enum_up, enum_dn);
    });
  }

  check_dimension_reasonable(size());
  check_dimension_works_with_blas_int_size(size());
}
XDIAG_CATCH

Electron::Electron(int64_t nsites) try
    : Electron(nsites, RepresentationSet{}) {}
XDIAG_CATCH

Electron::Electron(int64_t nsites, int64_t nup, int64_t ndn) try
    : Electron(nsites, RepresentationSet{Representation("nup", nup),
                                         Representation("ndn", ndn)}) {}
XDIAG_CATCH

Electron::Electron(int64_t nsites, Representation const &irrep) try
    : Electron(nsites, RepresentationSet{irrep}) {}
XDIAG_CATCH

Electron::Electron(int64_t nsites, int64_t nup, int64_t ndn,
                   Representation const &irrep) try
    : Electron(nsites, RepresentationSet{Representation("nup", nup),
                                         Representation("ndn", ndn), irrep}) {}
XDIAG_CATCH

int64_t Electron::nsites() const { return basis_->nsites(); }
int64_t Electron::index(ProductState const &pstate) const {
  return basis_->index(pstate);
}
int64_t Electron::dim() const { return size(); }
int64_t Electron::size() const { return basis_->size(); }
bool Electron::isreal() const { return irreps_.isreal(); }

ElectronIterator Electron::begin() const { return {this, 0}; }
ElectronIterator Electron::end() const { return {this, size()}; }

bool Electron::operator==(Electron const &rhs) const {
  return (irreps_ == rhs.irreps_) && (basis_ == rhs.basis_);
}

bool Electron::operator!=(Electron const &rhs) const {
  return !operator==(rhs);
}

RepresentationSet Electron::irreps() const { return irreps_; }
std::shared_ptr<basis::Basis> const &Electron::basis() const { return basis_; }

std::ostream &operator<<(std::ostream &out, Electron const &block) {
  print_block(out, block);
  return out;
}
std::string to_string(Electron const &block) {
  return utils::to_string_generic(block);
}

ElectronIterator::ElectronIterator(Electron const *block, int64_t idx)
    : idx_(idx) {
  if (idx_ < block->size()) {
    it_ = block->basis()->product_state_iterator();
  }
}

ElectronIterator &ElectronIterator::operator++() {
  it_->advance();
  ++idx_;
  return *this;
}

ProductState ElectronIterator::operator*() const {
  return it_->product_state();
}

bool ElectronIterator::operator==(ElectronIterator const &rhs) const {
  return idx_ == rhs.idx_;
}
bool ElectronIterator::operator!=(ElectronIterator const &rhs) const {
  return idx_ != rhs.idx_;
}

} // namespace xdiag
