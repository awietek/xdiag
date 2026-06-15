// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "spinhalf.hpp"

#include <regex>
#include <type_traits>

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_sublattice.hpp>
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
#include <xdiag/extern/fmt/color.hpp>
#endif

namespace xdiag {

namespace {

// Carries a type into a generic lambda without constructing a value of it.
template <typename T> struct type_tag {
  using type = T;
};

// Selects the (bit type, enumeration type) for a Spinhalf block from nsites and
// the optional magnetization nup, builds the enumeration, and forwards it to f
// as f(enumeration). The single place the nsites -> bit/enum ladder lives;
// BasisOnTheFly and BasisSymmetric both go through it.
template <typename F>
void dispatch_enumeration(int64_t nsites, std::optional<int64_t> nup, F &&f) {
  using namespace bits;
  using namespace combinatorics;
  if (!nup) { // full Hilbert space
    if (nsites <= 32) {
      f(Subsets<uint32_t>(nsites));
    } else if (nsites <= 64) {
      f(Subsets<uint64_t>(nsites));
    } else {
      XDIAG_THROW(
          "Unsupported nsites > 64 for non-Sz conserving Spinhalf block");
    }
  } else { // fixed magnetization; LinTable for fast lookup up to 42 sites
    int64_t n = nsites;
    int64_t k = *nup;
    if (nsites <= 32) {
      f(LinTable<uint32_t>(n, k));
    } else if (nsites <= 42) {
      f(LinTable<uint64_t>(n, k));
    } else if (nsites <= 64) {
      f(Combinations<uint64_t>(n, k));
    } else if (nsites <= 128) {
      f(Combinations<BitsetStatic2>(n, k));
    } else if (nsites <= 256) {
      f(Combinations<BitsetStatic4>(n, k));
    } else if (nsites <= 512) {
      f(Combinations<BitsetStatic8>(n, k));
    } else {
      f(Combinations<BitsetDynamic>(n, k));
    }
  }
}

// Selects the BasisSublattice type from a "<k>sublattice" backend (k in 1..5)
// and nsites, and forwards it as a type_tag to f.
template <int n_sublat = 1, typename F>
void dispatch_sublattice(std::string const &backend, int64_t nsites, F &&f) {
  using namespace basis;
  if constexpr (n_sublat <= 5) {
    if (backend != fmt::format("{}sublattice", n_sublat)) {
      dispatch_sublattice<n_sublat + 1>(backend, nsites, f);
    } else if (nsites <= 32) {
      f(type_tag<BasisSublattice<uint32_t, n_sublat>>{});
    } else if (nsites <= 64) {
      f(type_tag<BasisSublattice<uint64_t, n_sublat>>{});
    } else {
      XDIAG_THROW(
          "Unsupported nsites > 64 for Spinhalf block with sublattice coding");
    }
  } else {
    XDIAG_THROW(fmt::format("Unknown backend: \"{}\"", backend));
  }
}

} // namespace

Spinhalf::Spinhalf(int64_t nsites, RepresentationSet const &irreps,
                   std::string backend) try
    : irreps_(irreps) {
  using namespace basis;

  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  }

  std::optional<int64_t> nup = irreps.charge("nup");
  if (nup) {
    if (*nup < 0) {
      XDIAG_THROW("Invalid argument: nup < 0");
    } else if (*nup > nsites) {
      XDIAG_THROW("Invalid argument: nup > nsites");
    }
  }

  std::optional<PermutationGroup> group = irreps.group("SitePermutation");
  std::optional<Vector> characters = irreps.characters("SitePermutation");

  if (!group) { // no permutation symmetry -> BasisOnTheFly
    dispatch_enumeration(nsites, nup, [&](auto &&enumeration) {
      using enum_t = std::decay_t<decltype(enumeration)>;
      basis_ = std::make_shared<BasisOnTheFly<enum_t>>(enumeration);
    });
  } else if (backend == "auto") { // permutation symmetry -> BasisSymmetric
    dispatch_enumeration(nsites, nup, [&](auto &&enumeration) {
      using enum_t = std::decay_t<decltype(enumeration)>;
      basis_ = std::make_shared<BasisSymmetric<enum_t>>(enumeration, *group,
                                                        *characters);
    });
  } else { // "<k>sublattice" coding -> BasisSublattice
    dispatch_sublattice(backend, nsites, [&](auto tag) {
      using basis_t = typename decltype(tag)::type;
      if (nup) {
        basis_ = std::make_shared<basis_t>(*nup, *group, *characters);
      } else {
        basis_ = std::make_shared<basis_t>(*group, *characters);
      }
    });
  }

  check_dimension_reasonable(size());
  check_dimension_works_with_blas_int_size(size());
}
XDIAG_CATCH

Spinhalf::Spinhalf(int64_t nsites) try
    : Spinhalf(nsites, RepresentationSet{}) {}
XDIAG_CATCH

Spinhalf::Spinhalf(int64_t nsites, int64_t nup) try
    : Spinhalf(nsites, RepresentationSet{Representation("nup", nup)}) {}
XDIAG_CATCH

Spinhalf::Spinhalf(int64_t nsites, Representation const &irrep,
                   std::string backend) try
    : Spinhalf(nsites, RepresentationSet{irrep}, backend) {}
XDIAG_CATCH

Spinhalf::Spinhalf(int64_t nsites, int64_t nup, Representation const &irrep,
                   std::string backend) try
    : Spinhalf(nsites, RepresentationSet{Representation("nup", nup), irrep},
               backend) {}
XDIAG_CATCH

int64_t Spinhalf::nsites() const { return basis_->nsites(); }
int64_t Spinhalf::dim() const { return size(); }
int64_t Spinhalf::size() const { return basis_->size(); }
bool Spinhalf::isreal() const { return irreps_.isreal(); }

SpinhalfIterator Spinhalf::begin() const { return {this, 0}; }
SpinhalfIterator Spinhalf::end() const { return {this, size()}; }

bool Spinhalf::operator==(Spinhalf const &rhs) const {
  return (irreps_ == rhs.irreps_) && (basis_ == rhs.basis_);
}

bool Spinhalf::operator!=(Spinhalf const &rhs) const {
  return !operator==(rhs);
}

RepresentationSet Spinhalf::irreps() const { return irreps_; }
std::shared_ptr<basis::Basis> const &Spinhalf::basis() const { return basis_; }

std::ostream &operator<<(std::ostream &out, Spinhalf const &block) {
  print_block(out, block);
  return out;
}
std::string to_string(Spinhalf const &block) {
  return to_string_generic(block);
}

SpinhalfIterator::SpinhalfIterator(Spinhalf const *block, int64_t idx)
    : idx_(idx) {
  if (idx_ < block->size()) {
    it_ = block->basis()->product_state_iterator();
  }
}

SpinhalfIterator &SpinhalfIterator::operator++() {
  it_->advance();
  ++idx_;
  return *this;
}

ProductState SpinhalfIterator::operator*() const {
  return it_->product_state();
}

bool SpinhalfIterator::operator==(SpinhalfIterator const &rhs) const {
  return idx_ == rhs.idx_;
}
bool SpinhalfIterator::operator!=(SpinhalfIterator const &rhs) const {
  return idx_ != rhs.idx_;
}

} // namespace xdiag
