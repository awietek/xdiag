// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fermion.hpp"

#include <regex>

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/blocks/print_block.hpp>
#include <xdiag/random/hash.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>

#ifndef XDIAG_DISABLE_COLOR
#include <extern/fmt/color.hpp>
#endif

namespace xdiag {

template <typename F>
static void dispatch_enumeration(int64_t nsites, std::optional<int64_t> number,
                                 F &&f) {
  using namespace bits;
  using namespace combinatorics;
  if (number) { // fixed magnetization; LinTable for fast lookup up to 42 sites
    int64_t n = nsites;
    int64_t k = *number;
    if (nsites <= 32) {
      f(LinTable<uint32_t>(n, k));
    } else if (nsites <= 42) {
      f(LinTable<uint64_t>(n, k));
    } else if (nsites <= 64) {
      f(Combinations<uint64_t>(n, k));
    } else if (nsites <= 128) {
      f(Combinations<BitsetStatic2>(n, k));
    } else {
      f(Combinations<BitsetDynamic>(n, k));
    }
  } else { // full Hilbert space
    if (nsites <= 32) {
      f(Subsets<uint32_t>(nsites));
    } else if (nsites <= 64) {
      f(Subsets<uint64_t>(nsites));
    } else {
      XDIAG_THROW(
          "Unsupported nsites > 64 for non-number conserving Fermion block");
    }
  }
}

Fermion::Fermion(int64_t nsites, RepresentationSet const &irreps) try
    : irreps_(irreps) {
  using namespace basis;

  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  }

  std::optional<int64_t> number = irreps.charge("number");
  if (number) {
    if (*number < 0) {
      XDIAG_THROW("Invalid argument: number < 0");
    } else if (*number > nsites) {
      XDIAG_THROW("Invalid argument: number > nsites");
    }
  }

  std::optional<PermutationGroup> group = irreps.group("SitePermutation");
  std::optional<Vector> characters = irreps.characters("SitePermutation");

  if (group) { // permutation symmetry -> BasisSymmetric (fermionic)
    dispatch_enumeration(nsites, number, [&](auto &&enumeration) {
      using enum_t = std::decay_t<decltype(enumeration)>;
      basis_ = std::make_shared<BasisSymmetric<enum_t>>(
          enumeration, *group, *characters, /*fermionic=*/true);
    });
  } else { // no permutation symmetry -> BasisOnTheFly
    dispatch_enumeration(nsites, number, [&](auto &&enumeration) {
      using enum_t = std::decay_t<decltype(enumeration)>;
      basis_ = std::make_shared<BasisOnTheFly<enum_t>>(enumeration);
    });
  }

  check_dimension_reasonable(size());
  check_dimension_works_with_blas_int_size(size());
}
XDIAG_CATCH

Fermion::Fermion(int64_t nsites) try : Fermion(nsites, RepresentationSet{}) {}
XDIAG_CATCH

Fermion::Fermion(int64_t nsites, int64_t number) try
    : Fermion(nsites, RepresentationSet{Representation("number", number)}) {}
XDIAG_CATCH

Fermion::Fermion(int64_t nsites, Representation const &irrep) try
    : Fermion(nsites, RepresentationSet{irrep}) {}
XDIAG_CATCH

Fermion::Fermion(int64_t nsites, int64_t number,
                 Representation const &irrep) try
    : Fermion(nsites,
              RepresentationSet{Representation("number", number), irrep}) {}
XDIAG_CATCH

int64_t Fermion::nsites() const { return basis_->nsites(); }
int64_t Fermion::dim() const { return size(); }
int64_t Fermion::size() const { return basis_->size(); }
bool Fermion::isreal() const { return irreps_.isreal(); }

FermionIterator Fermion::begin() const { return {this, 0}; }
FermionIterator Fermion::end() const { return {this, size()}; }

bool Fermion::operator==(Fermion const &rhs) const {
  return (irreps_ == rhs.irreps_) && (basis_ == rhs.basis_);
}

bool Fermion::operator!=(Fermion const &rhs) const { return !operator==(rhs); }

RepresentationSet Fermion::irreps() const { return irreps_; }
std::shared_ptr<basis::Basis> const &Fermion::basis() const { return basis_; }

std::ostream &operator<<(std::ostream &out, Fermion const &block) {
  print_block(out, block);
  return out;
}
std::string to_string(Fermion const &block) { return to_string_generic(block); }

FermionIterator::FermionIterator(Fermion const *block, int64_t idx)
    : idx_(idx) {
  if (idx_ < block->size()) {
    it_ = block->basis()->product_state_iterator();
  }
}

FermionIterator &FermionIterator::operator++() {
  it_->advance();
  ++idx_;
  return *this;
}

ProductState FermionIterator::operator*() const { return it_->product_state(); }

bool FermionIterator::operator==(FermionIterator const &rhs) const {
  return idx_ == rhs.idx_;
}
bool FermionIterator::operator!=(FermionIterator const &rhs) const {
  return idx_ != rhs.idx_;
}

} // namespace xdiag
