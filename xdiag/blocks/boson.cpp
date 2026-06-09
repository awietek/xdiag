// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "boson.hpp"

#include <memory>
#include <optional>
#include <type_traits>
#include <utility>

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/bounded_multisets/bounded_multisets.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>
#include <xdiag/combinatorics/bounded_partitions/schaefer_table.hpp>
#include <xdiag/math/log2.hpp>
#include <xdiag/random/hash.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>

#ifndef XDIAG_DISABLE_COLOR
#include <xdiag/extern/fmt/color.hpp>
#endif

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

// Dispatches the runtime number of bits per site (nlocalbits in 1..8) to the
// compile-time template parameter N.
template <typename bit_t, typename F>
void dispatch_nlocalbits(int64_t nlocalbits, int64_t nsites, int64_t d,
                         std::optional<int64_t> number, F &&f) try {
  switch (nlocalbits) {
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
  case 5:
    dispatch_boson_enumeration<bit_t, 5>(nsites, d, number, std::forward<F>(f));
    break;
  case 6:
    dispatch_boson_enumeration<bit_t, 6>(nsites, d, number, std::forward<F>(f));
    break;
  case 7:
    dispatch_boson_enumeration<bit_t, 7>(nsites, d, number, std::forward<F>(f));
    break;
  case 8:
    dispatch_boson_enumeration<bit_t, 8>(nsites, d, number, std::forward<F>(f));
    break;
  default:
    XDIAG_THROW(fmt::format("Unsupported number of bits per site: {}. "
                            "nlocalbits must lie in 1..8.",
                            nlocalbits));
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
  if (nsites * nlocalbits <= 64) {
    dispatch_nlocalbits<uint64_t>(nlocalbits, nsites, d, number,
                                  std::forward<F>(f));
  } else {
    dispatch_nlocalbits<BitsetDynamic>(nlocalbits, nsites, d, number,
                                       std::forward<F>(f));
  }
}
XDIAG_CATCH

Boson::Boson(int64_t nsites, int64_t d, RepresentationSet const &irreps) try
    : nsites_(nsites), d_(d), irreps_(irreps) {
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

int64_t Boson::nsites() const { return nsites_; }
int64_t Boson::d() const { return d_; }
int64_t Boson::dim() const { return size(); }
int64_t Boson::size() const { return basis_->size(); }
bool Boson::isreal() const { return irreps_.isreal(); }

BosonIterator Boson::begin() const { return {this, 0}; }
BosonIterator Boson::end() const { return {this, size()}; }

bool Boson::operator==(Boson const &rhs) const {
  return (nsites_ == rhs.nsites_) && (d_ == rhs.d_) && (irreps_ == rhs.irreps_);
}

bool Boson::operator!=(Boson const &rhs) const { return !operator==(rhs); }

RepresentationSet Boson::irreps() const { return irreps_; }
std::shared_ptr<basis::Basis> const &Boson::basis() const { return basis_; }

int64_t nsites(Boson const &block) { return block.nsites(); }
int64_t dim(Boson const &block) { return block.dim(); }
int64_t size(Boson const &block) { return block.size(); }
bool isreal(Boson const &block) { return block.isreal(); }

std::ostream &operator<<(std::ostream &out, Boson const &block) {
  out << fmt::format(fg(fmt::color::steel_blue) | fmt::emphasis::bold,
                     "Boson\n");
  out << "| nsites   : " << block.nsites() << "\n";
  out << "| d        : " << block.d() << "\n";
  std::optional<int64_t> number = block.irreps().charge("number");
  if (number) {
    out << "| number   : " << *number << "\n";
  } else {
    out << "| number   : not conserved\n";
  }
  auto group = block.irreps().group("SitePermutation");
  if (group) {
    out << "| permutation symmetries used\n";
    out << fmt::format(
        "| irrep ID : {0:x}\n",
        random::hash(Representation(
            *group, *block.irreps().characters("SitePermutation"))));
  }
  out << fmt::format("| ID       : {0:x}\n", random::hash(block));
  out << "| dimension: " << fmt::format_de("{:L}", block.size()) << "\n";
  return out;
}
std::string to_string(Boson const &block) { return to_string_generic(block); }

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
