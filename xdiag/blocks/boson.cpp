// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "boson.hpp"

#include <memory>
#include <optional>
#include <regex>
#include <type_traits>
#include <utility>

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/blocks/print_block.hpp>
#include <xdiag/combinatorics/bounded_multisets/bounded_multisets.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>
#include <xdiag/combinatorics/bounded_partitions/schaefer_table.hpp>
#include <xdiag/math/log2.hpp>
#include <xdiag/random/hash.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>

namespace xdiag {

// Builds the boson enumeration for a fixed bit storage type (bit_t) and a fixed
// number of bits per site (N), then forwards it to f as f(enumeration). The
// enumerator is chosen at runtime:
//   number present -> SchaeferTable while the state fits a fast lookup table
//                     (<= 42 bits), otherwise BoundedPartitions (on-the-fly
//                     ranking). Both enumerate states at fixed particle number.
//   number absent  -> BoundedMultisets (all fillings up to d-1).
template <typename bit_t, int N, typename F>
void dispatch_boson_enumeration(int64_t nsites, int64_t d,
                                std::optional<int64_t> number, F &&f) try {
  using namespace bits;
  using namespace combinatorics;
  using bitarray_t = BitArray<bit_t, N>;

  if (number) {
    if constexpr (std::is_same_v<bit_t, uint64_t>) {
      if (nsites * N <= 42) {
        f(SchaeferTable<bitarray_t>(nsites, *number, d));
      } else {
        f(BoundedPartitions<bitarray_t>(nsites, *number, d));
      }
    } else {
      f(BoundedPartitions<bitarray_t>(nsites, *number, d));
    }
  } else {
    f(BoundedMultisets<bitarray_t>(nsites, d));
  }
}
XDIAG_CATCH

// Rounds the required bits-per-site up to the nearest compiled BitArray width.
// Only widths {1,2,3,4,8} are instantiated to keep the library small; widths
// 5,6,7 are stored in 8-bit fields instead. This is functionally identical (the
// enumeration is driven by the true local dimension d, N is only the storage
// field width), it just uses slightly wider storage per site.
constexpr int64_t promote_nlocalbits(int64_t nlocalbits) {
  return (nlocalbits <= 4) ? nlocalbits : 8;
}

// Dispatches the (already promoted) storage width N in {1,2,3,4,8} to the
// compile-time template parameter.
template <typename bit_t, typename F>
void dispatch_nlocalbits(int64_t nstore, int64_t nsites, int64_t d,
                         std::optional<int64_t> number, F &&f) try {
  switch (nstore) {
  case 1:
    dispatch_boson_enumeration<bit_t, 1>(nsites, d, number, std::forward<F>(f));
    break;
  case 2:
    dispatch_boson_enumeration<bit_t, 2>(nsites, d, number, std::forward<F>(f));
    break;
  case 3:
    dispatch_boson_enumeration<bit_t, 3>(nsites, d, number, std::forward<F>(f));
    break;
  case 4:
    dispatch_boson_enumeration<bit_t, 4>(nsites, d, number, std::forward<F>(f));
    break;
  case 8:
    dispatch_boson_enumeration<bit_t, 8>(nsites, d, number, std::forward<F>(f));
    break;
  default:
    XDIAG_THROW(fmt::format("Unsupported number of bits per site: {}. "
                            "nlocalbits must lie in 1..8.",
                            nstore));
  }
}
XDIAG_CATCH

// Selects the bit storage type (a single 64-bit word while the whole state fits
// in 64 bits, otherwise a dynamically sized bitset) and forwards the resulting
// boson enumeration to f. The single place the boson nsites/d -> bit/enum
// ladder lives; BasisOnTheFly and BasisSymmetric both go through it.
template <typename F>
void dispatch_boson_enumeration(int64_t nlocalbits, int64_t nsites, int64_t d,
                                std::optional<int64_t> number, F &&f) try {
  using namespace bits;
  // Promote to the actual storage width first, so the 64-bit fit test below
  // reflects how many bits are really used (nsites * nstore), not the ideal
  // width. Getting this wrong would pick a single-word backend for a state that
  // no longer fits in 64 bits.
  int64_t nstore = promote_nlocalbits(nlocalbits);
  if (nsites * nstore <= 64) {
    dispatch_nlocalbits<uint64_t>(nstore, nsites, d, number,
                                  std::forward<F>(f));
  } else {
    dispatch_nlocalbits<BitsetDynamic>(nstore, nsites, d, number,
                                       std::forward<F>(f));
  }
}
XDIAG_CATCH

Boson::Boson(int64_t nsites, int64_t d, RepresentationSet const &irreps) try
    : d_(d), irreps_(irreps) {
  using namespace basis;

  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  }

  if (d < 2) {
    XDIAG_THROW("Unsupported local dimension: d < 2");
  }

  int64_t nlocalbits = math::ceillog2(d);
  if (nlocalbits > 8) {
    XDIAG_THROW("Unsupported local dimension: d > 256");
  }

  std::optional<int64_t> number = irreps.charge("number");
  if (number) {
    if (*number < 0) {
      XDIAG_THROW("Invalid argument: number < 0");
    } else if (*number > nsites * (d - 1)) {
      XDIAG_THROW("Invalid argument: number > nsites * (d-1)");
    }
  }

  std::optional<PermutationGroup> group = irreps.group("SitePermutation");
  std::optional<Vector> characters = irreps.characters("SitePermutation");

  if (!group) { // no permutation symmetry -> BasisOnTheFly
    dispatch_boson_enumeration(
        nlocalbits, nsites, d, number, [&](auto &&enumeration) {
          using enum_t = std::decay_t<decltype(enumeration)>;
          basis_ = std::make_shared<BasisOnTheFly<enum_t>>(enumeration);
        });
  } else { // permutation symmetry -> BasisSymmetric
    dispatch_boson_enumeration(
        nlocalbits, nsites, d, number, [&](auto &&enumeration) {
          using enum_t = std::decay_t<decltype(enumeration)>;
          basis_ = std::make_shared<BasisSymmetric<enum_t>>(enumeration, *group,
                                                            *characters);
        });
  }

  check_dimension_reasonable(size());
  check_dimension_works_with_blas_int_size(size());
}
XDIAG_CATCH

Boson::Boson(int64_t nsites, int64_t d) try
    : Boson(nsites, d, RepresentationSet{}) {}
XDIAG_CATCH

Boson::Boson(int64_t nsites, int64_t d, int64_t number) try
    : Boson(nsites, d, RepresentationSet{Representation("number", number)}) {}
XDIAG_CATCH

Boson::Boson(int64_t nsites, int64_t d, Representation const &irrep) try
    : Boson(nsites, d, RepresentationSet{irrep}) {}
XDIAG_CATCH

Boson::Boson(int64_t nsites, int64_t d, int64_t number,
             Representation const &irrep) try
    : Boson(nsites, d,
            RepresentationSet{Representation("number", number), irrep}) {}
XDIAG_CATCH

int64_t Boson::nsites() const { return basis_->nsites(); }
int64_t Boson::index(ProductState const &pstate) const {
  return basis_->index(pstate);
}
int64_t Boson::d() const { return d_; }
int64_t Boson::dim() const { return size(); }
int64_t Boson::size() const { return basis_->size(); }
bool Boson::isreal() const { return irreps_.isreal(); }

BosonIterator Boson::begin() const { return {this, 0}; }
BosonIterator Boson::end() const { return {this, size()}; }

bool Boson::operator==(Boson const &rhs) const {
  return (d_ == rhs.d_) && (irreps_ == rhs.irreps_) && (basis_ == rhs.basis_);
}

bool Boson::operator!=(Boson const &rhs) const { return !operator==(rhs); }

RepresentationSet Boson::irreps() const { return irreps_; }
std::shared_ptr<basis::Basis> const &Boson::basis() const { return basis_; }

std::ostream &operator<<(std::ostream &out, Boson const &block) {
  print_block(out, block);
  return out;
}
std::string to_string(Boson const &block) {
  return utils::to_string_generic(block);
}

BosonIterator::BosonIterator(Boson const *block, int64_t idx) : idx_(idx) {
  if (idx_ < block->size()) {
    it_ = block->basis()->product_state_iterator();
  }
}

BosonIterator &BosonIterator::operator++() {
  it_->advance();
  ++idx_;
  return *this;
}

ProductState BosonIterator::operator*() const { return it_->product_state(); }

bool BosonIterator::operator==(BosonIterator const &rhs) const {
  return idx_ == rhs.idx_;
}
bool BosonIterator::operator!=(BosonIterator const &rhs) const {
  return idx_ != rhs.idx_;
}

} // namespace xdiag
