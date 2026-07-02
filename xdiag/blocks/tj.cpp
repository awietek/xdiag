// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "tj.hpp"

#include <optional>

#include <xdiag/basis/basis_tj.hpp>
#include <xdiag/basis/basis_tj_symmetric.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/blocks/print_block.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/math/vector.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>

namespace xdiag {

tJ::tJ(int64_t nsites, RepresentationSet const &irreps) try : irreps_(irreps) {
  using namespace basis;
  using namespace combinatorics;
  using namespace bits;

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
  if (nup && ndn && ((*nup + *ndn) > nsites)) {
    XDIAG_THROW("Invalid argument: nup + ndn > nsites (no double occupancy)");
  }

  std::optional<PermutationGroup> group = irreps.group("SitePermutation");
  std::optional<Vector> characters = irreps.characters("SitePermutation");

  if (group) {
    // permutation symmetry -> BasistJSymmetric. Unlike the non-symmetric block
    // (compressed dn on nsites-nup sites), the symmetric basis stores FULL-nsites
    // dn, so enum_dn is built on (nsites, ndn) -- not (nsites-nup, ndn).
    auto build = [&](auto &&enum_up, auto &&enum_dn) {
      using enum_t = std::decay_t<decltype(enum_up)>;
      basis_ = std::make_shared<BasistJSymmetric<enum_t>>(enum_up, enum_dn,
                                                          *group, *characters);
    };
    if (nup && ndn) {
      int64_t n = nsites, ku = *nup, kd = *ndn;
      if (nsites <= 32) {
        build(LinTable<uint32_t>(n, ku), LinTable<uint32_t>(n, kd));
      } else if (nsites <= 42) {
        build(LinTable<uint64_t>(n, ku), LinTable<uint64_t>(n, kd));
      } else if (nsites <= 64) {
        build(Combinations<uint64_t>(n, ku), Combinations<uint64_t>(n, kd));
      } else {
        build(Combinations<BitsetDynamic>(n, ku),
              Combinations<BitsetDynamic>(n, kd));
      }
    } else if (!nup && !ndn) {
      if (nsites <= 32) {
        build(Subsets<uint32_t>(nsites), Subsets<uint32_t>(nsites));
      } else if (nsites <= 64) {
        build(Subsets<uint64_t>(nsites), Subsets<uint64_t>(nsites));
      } else {
        XDIAG_THROW("Unsupported nsites > 64 for non-number conserving "
                    "symmetric tJ block");
      }
    } else {
      XDIAG_THROW("tJ block requires both nup and ndn to be set, or neither");
    }
  } else if (nup && ndn) { // number-conserving: ups on (nsites,nup), compressed
                           // dn on (nsites-nup, ndn)
    int64_t n = nsites, ku = *nup, kd = *ndn, nc = nsites - *nup;
    auto build = [&](auto &&enum_up, auto &&enum_dncs) {
      using enum_t = std::decay_t<decltype(enum_up)>;
      basis_ = std::make_shared<BasistJ<enum_t>>(enum_up, enum_dncs);
    };
    if (nsites <= 32) {
      build(LinTable<uint32_t>(n, ku), LinTable<uint32_t>(nc, kd));
    } else if (nsites <= 42) {
      build(LinTable<uint64_t>(n, ku), LinTable<uint64_t>(nc, kd));
    } else if (nsites <= 64) {
      build(Combinations<uint64_t>(n, ku), Combinations<uint64_t>(nc, kd));
    } else {
      build(Combinations<BitsetDynamic>(n, ku),
            Combinations<BitsetDynamic>(nc, kd));
    }
  } else if (!nup && !ndn) { // no number conservation (small nsites only)
    auto build = [&](auto &&enum_up) {
      using enum_t = std::decay_t<decltype(enum_up)>;
      basis_ = std::make_shared<BasistJ<enum_t>>(enum_up);
    };
    if (nsites <= 32) {
      build(Subsets<uint32_t>(nsites));
    } else if (nsites <= 64) {
      build(Subsets<uint64_t>(nsites));
    } else {
      XDIAG_THROW("Unsupported nsites > 64 for non-number conserving tJ block");
    }
  } else {
    XDIAG_THROW("tJ block requires both nup and ndn to be set, or neither");
  }

  check_dimension_reasonable(size());
  check_dimension_works_with_blas_int_size(size());
}
XDIAG_CATCH

tJ::tJ(int64_t nsites) try : tJ(nsites, RepresentationSet{}) {}
XDIAG_CATCH

tJ::tJ(int64_t nsites, int64_t nup, int64_t ndn) try
    : tJ(nsites, RepresentationSet{Representation("nup", nup),
                                   Representation("ndn", ndn)}) {}
XDIAG_CATCH

tJ::tJ(int64_t nsites, Representation const &irrep) try
    : tJ(nsites, RepresentationSet{irrep}) {}
XDIAG_CATCH

tJ::tJ(int64_t nsites, int64_t nup, int64_t ndn, Representation const &irrep) try
    : tJ(nsites, RepresentationSet{Representation("nup", nup),
                                   Representation("ndn", ndn), irrep}) {}
XDIAG_CATCH

int64_t tJ::nsites() const { return basis_->nsites(); }
int64_t tJ::index(ProductState const &pstate) const {
  return basis_->index(pstate);
}
int64_t tJ::dim() const { return size(); }
int64_t tJ::size() const { return basis_->size(); }
bool tJ::isreal() const { return irreps_.isreal(); }

tJIterator tJ::begin() const { return {this, 0}; }
tJIterator tJ::end() const { return {this, size()}; }

bool tJ::operator==(tJ const &rhs) const {
  return (irreps_ == rhs.irreps_) && (basis_ == rhs.basis_);
}
bool tJ::operator!=(tJ const &rhs) const { return !operator==(rhs); }

RepresentationSet tJ::irreps() const { return irreps_; }
std::shared_ptr<basis::Basis> const &tJ::basis() const { return basis_; }

std::ostream &operator<<(std::ostream &out, tJ const &block) {
  print_block(out, block);
  return out;
}
std::string to_string(tJ const &block) { return to_string_generic(block); }

tJIterator::tJIterator(tJ const *block, int64_t idx) : idx_(idx) {
  if (idx_ < block->size()) {
    it_ = block->basis()->product_state_iterator();
  }
}

tJIterator &tJIterator::operator++() {
  it_->advance();
  ++idx_;
  return *this;
}

ProductState tJIterator::operator*() const { return it_->product_state(); }

bool tJIterator::operator==(tJIterator const &rhs) const {
  return idx_ == rhs.idx_;
}
bool tJIterator::operator!=(tJIterator const &rhs) const {
  return idx_ != rhs.idx_;
}

} // namespace xdiag
