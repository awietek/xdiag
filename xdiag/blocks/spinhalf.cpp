// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "spinhalf.hpp"

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/math/binomial.hpp>
#include <xdiag/math/ipow.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>

namespace xdiag {

Spinhalf::Spinhalf(int64_t nsites) try
    : nsites_(nsites), nup_(std::nullopt), irrep_(std::nullopt) {
  using namespace bits;
  using namespace combinatorics;
  using namespace basis;

  // Safety check
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  }

  size_ = math::ipow(2, nsites);
  check_dimension_reasonable(size_);
  check_dimension_works_with_blas_int_size(size_);

  // Choose basis implementation
  if (nsites <= 32) {
    using bit_t = uint32_t;
    using enum_t = Subsets<bit_t>;
    using basis_t = BasisOnTheFly<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites));
  } else if (nsites <= 64) {
    using bit_t = uint64_t;
    using enum_t = Subsets<bit_t>;
    using basis_t = BasisOnTheFly<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites));
  } else {
    XDIAG_THROW("Invalid nsites > 64 for non-Sz conserving Spinhalf block");
  }
}
XDIAG_CATCH

Spinhalf::Spinhalf(int64_t nsites, int64_t nup) try
    : nsites_(nsites), nup_(nup), irrep_(std::nullopt) {
  using namespace bits;
  using namespace combinatorics;
  using namespace basis;

  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  } else if (nup < 0) {
    XDIAG_THROW("Invalid argument: nup < 0");
  } else if (nup > nsites) {
    XDIAG_THROW("Invalid argument: nup > nsites");
  }

  size_ = math::binomial(nsites, nup);
  check_dimension_reasonable(size_);
  check_dimension_works_with_blas_int_size(size_);

  // For nsites <= 42 choose a LinTable for fast lookups
  if (nsites <= 32) {
    using bit_t = uint32_t;
    using enum_t = LinTable<bit_t>;
    using basis_t = BasisOnTheFly<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup));
  } else if (nsites <= 42) {
    using bit_t = uint64_t;
    using enum_t = LinTable<bit_t>;
    using basis_t = BasisOnTheFly<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup));
  } else if (nsites <= 64) {
    using bit_t = uint64_t;
    using enum_t = Combinations<bit_t>;
    using basis_t = BasisOnTheFly<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup));
  } else if (nsites <= 128) {
    using bit_t = Bitset<uint64_t, 2>;
    using enum_t = Combinations<bit_t>;
    using basis_t = BasisOnTheFly<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup));
  } else if (nsites <= 256) {
    using bit_t = Bitset<uint64_t, 4>;
    using enum_t = Combinations<bit_t>;
    using basis_t = BasisOnTheFly<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup));
  } else if (nsites <= 512) {
    using bit_t = Bitset<uint64_t, 8>;
    using enum_t = Combinations<bit_t>;
    using basis_t = BasisOnTheFly<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup));
  } else { // use dynamically sized Bitset
    using bit_t = Bitset<uint64_t, 0>;
    using enum_t = Combinations<bit_t>;
    using basis_t = BasisOnTheFly<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup));
  }
}
XDIAG_CATCH

Spinhalf::Spinhalf(int64_t nsites, Representation const &irrep) try
    : nsites_(nsites), nup_(std::nullopt), irrep_(irrep) {
  using namespace bits;
  using namespace combinatorics;
  using namespace basis;

  // Safety check
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  }

  // Choose basis implementation
  if (nsites <= 32) {
    using bit_t = uint32_t;
    using enum_t = Subsets<bit_t>;
    using basis_t = BasisSymmetric<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites), irrep);
  } else if (nsites <= 64) {
    using bit_t = uint64_t;
    using enum_t = Subsets<bit_t>;
    using basis_t = BasisSymmetric<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites), irrep);
  } else {
    XDIAG_THROW("Invalid nsites > 64 for non-Sz conserving Spinhalf block");
  }

  size_ = basis_->size();
  check_dimension_reasonable(size_);
  check_dimension_works_with_blas_int_size(size_);
}
XDIAG_CATCH

Spinhalf::Spinhalf(int64_t nsites, int64_t nup, Representation const &irrep) try
    : nsites_(nsites), nup_(nup), irrep_(irrep) {
  using namespace bits;
  using namespace combinatorics;
  using namespace basis;

  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  } else if (nup < 0) {
    XDIAG_THROW("Invalid argument: nup < 0");
  } else if (nup > nsites) {
    XDIAG_THROW("Invalid argument: nup > nsites");
  }

  // For nsites <= 42 choose a LinTable for fast lookups
  if (nsites <= 32) {
    using bit_t = uint32_t;
    using enum_t = LinTable<bit_t>;
    using basis_t = BasisSymmetric<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup), irrep);
  } else if (nsites <= 42) {
    using bit_t = uint64_t;
    using enum_t = LinTable<bit_t>;
    using basis_t = BasisSymmetric<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup), irrep);
  } else if (nsites <= 64) {
    using bit_t = uint64_t;
    using enum_t = Combinations<bit_t>;
    using basis_t = BasisSymmetric<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup), irrep);
  } else if (nsites <= 128) {
    using bit_t = Bitset<uint64_t, 2>;
    using enum_t = Combinations<bit_t>;
    using basis_t = BasisSymmetric<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup), irrep);
  } else if (nsites <= 256) {
    using bit_t = Bitset<uint64_t, 4>;
    using enum_t = Combinations<bit_t>;
    using basis_t = BasisSymmetric<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup), irrep);
  } else if (nsites <= 512) {
    using bit_t = Bitset<uint64_t, 8>;
    using enum_t = Combinations<bit_t>;
    using basis_t = BasisSymmetric<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup), irrep);
  } else { // use dynamically sized Bitset
    using bit_t = Bitset<uint64_t, 0>;
    using enum_t = Combinations<bit_t>;
    using basis_t = BasisSymmetric<enum_t>;
    basis_ = std::make_shared<basis_t>(enum_t(nsites, nup), irrep);
  }
  size_ = basis_->size();
  check_dimension_reasonable(size_);
  check_dimension_works_with_blas_int_size(size_);
}
XDIAG_CATCH

int64_t Spinhalf::dim() const { return size_; }
int64_t Spinhalf::size() const { return size_; }
bool Spinhalf::isreal() const { return irrep_ ? irrep_->isreal() : true; }

SpinhalfIterator Spinhalf::begin() const { return {this, 0}; }
SpinhalfIterator Spinhalf::end() const { return {this, size_}; }

bool Spinhalf::operator==(Spinhalf const &rhs) const {
  return (nsites_ == rhs.nsites_) && (nup_ == rhs.nup_) &&
         (size_ == rhs.size_) && (basis_ == rhs.basis_);
}

bool Spinhalf::operator!=(Spinhalf const &rhs) const {
  return !operator==(rhs);
}

int64_t Spinhalf::nsites() const { return nsites_; }
std::optional<int64_t> Spinhalf::nup() const { return nup_; }
std::shared_ptr<basis::Basis> const &Spinhalf::basis() const { return basis_; }

int64_t nsites(Spinhalf const &block) { return block.nsites(); }
int64_t dim(Spinhalf const &block) { return block.dim(); }
int64_t size(Spinhalf const &block) { return block.size(); }
bool isreal(Spinhalf const &block) { return block.isreal(); }

std::ostream &operator<<(std::ostream &out, Spinhalf const &block) {
  out << "Spinhalf:\n";
  out << "  nsites   : " << block.nsites() << "\n";
  if (block.nup()) {
    out << "  nup      : " << *block.nup() << "\n";
  } else {
    out << "  nup      : not conserved\n";
  }

  out << "  dimension: " << fmt::format_de("{:L}", block.size()) << "\n";
  return out;
}
std::string to_string(Spinhalf const &block) {
  return to_string_generic(block);
}

SpinhalfIterator::SpinhalfIterator(Spinhalf const *block, int64_t idx)
    : block_(block), idx_(idx) {}

SpinhalfIterator &SpinhalfIterator::operator++() {
  ++idx_;
  return *this;
}

ProductState SpinhalfIterator::operator*() const {
  return block_->basis()->product_state(idx_, {"Dn", "Up"});
}

bool SpinhalfIterator::operator==(SpinhalfIterator const &rhs) const {
  return idx_ == rhs.idx_;
}
bool SpinhalfIterator::operator!=(SpinhalfIterator const &rhs) const {
  return idx_ != rhs.idx_;
}

} // namespace xdiag
