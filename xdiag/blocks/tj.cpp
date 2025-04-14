// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "tj.hpp"

#include <xdiag/random/hash.hpp>

namespace xdiag {

using namespace basis;

tJ::tJ(int64_t nsites, int64_t nup, int64_t ndn, std::string backend) try
    : nsites_(nsites), backend_(backend), nup_(nup), ndn_(ndn),
      irrep_(std::nullopt) {
  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  } else if ((nup < 0) || (ndn < 0)) {
    XDIAG_THROW("Invalid argument: (nup < 0) or (ndn < 0)");
  } else if ((nup + ndn) > nsites) {
    XDIAG_THROW("Invalid argument: nup + ndn > nsites");
  }

  // Choose basis implementation
  if (backend == "auto") {
    if (nsites < 32) {
      basis_ =
          std::make_shared<basis_t>(tj::BasisNp<uint32_t>(nsites, nup, ndn));
    } else if (nsites < 64) {
      basis_ =
          std::make_shared<basis_t>(tj::BasisNp<uint64_t>(nsites, nup, ndn));
    } else {
      XDIAG_THROW("blocks with more than 64 sites currently not implemented");
    }
  } else if (backend == "32bit") {
    basis_ = std::make_shared<basis_t>(tj::BasisNp<uint32_t>(nsites, nup, ndn));
  } else if (backend == "64bit") {
    basis_ = std::make_shared<basis_t>(tj::BasisNp<uint64_t>(nsites, nup, ndn));
  } else {
    XDIAG_THROW(fmt::format("Unknown backend: \"{}\"", backend));
  }
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

tJ::tJ(int64_t nsites, int64_t nup, int64_t ndn, Representation const &irrep,
       std::string backend) try
    : nsites_(nsites), backend_(backend), nup_(nup), ndn_(ndn), irrep_(irrep) {
  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  } else if ((nup < 0) || (ndn < 0)) {
    XDIAG_THROW("Invalid argument: (nup < 0) or (ndn < 0)");
  } else if ((nup + ndn) > nsites) {
    XDIAG_THROW("Invalid argument: nup + ndn > nsites");
  } else if (nsites != irrep.group().nsites()) {
    XDIAG_THROW("nsites does not match the nsites in PermutationGroup");
  }

  // Choose basis implementation
  if (backend == "auto") {
    if (nsites < 32) {
      if (irrep.isreal()) {
        auto characters = irrep.characters().as<arma::vec>();
        basis_ = std::make_shared<basis_t>(tj::BasisSymmetricNp<uint32_t>(
            nsites, nup, ndn, irrep.group(), characters));
      } else {
        auto characters = irrep.characters().as<arma::cx_vec>();
        basis_ = std::make_shared<basis_t>(tj::BasisSymmetricNp<uint32_t>(
            nsites, nup, ndn, irrep.group(), characters));
      }

    } else if (nsites < 64) {
      if (irrep.isreal()) {
        auto characters = irrep.characters().as<arma::vec>();
        basis_ = std::make_shared<basis_t>(tj::BasisSymmetricNp<uint64_t>(
            nsites, nup, ndn, irrep.group(), characters));
      } else {
        auto characters = irrep.characters().as<arma::cx_vec>();
        basis_ = std::make_shared<basis_t>(tj::BasisSymmetricNp<uint64_t>(
            nsites, nup, ndn, irrep.group(), characters));
      }
    } else {
      XDIAG_THROW("blocks with more than 64 sites currently not implemented");
    }
  } else if (backend == "32bit") {
    if (irrep.isreal()) {
      auto characters = irrep.characters().as<arma::vec>();
      basis_ = std::make_shared<basis_t>(tj::BasisSymmetricNp<uint32_t>(
          nsites, nup, ndn, irrep.group(), characters));
    } else {
      auto characters = irrep.characters().as<arma::cx_vec>();
      basis_ = std::make_shared<basis_t>(tj::BasisSymmetricNp<uint32_t>(
          nsites, nup, ndn, irrep.group(), characters));
    }
  } else if (backend == "64bit") {
    if (irrep.isreal()) {
      auto characters = irrep.characters().as<arma::vec>();
      basis_ = std::make_shared<basis_t>(tj::BasisSymmetricNp<uint64_t>(
          nsites, nup, ndn, irrep.group(), characters));
    } else {
      auto characters = irrep.characters().as<arma::cx_vec>();
      basis_ = std::make_shared<basis_t>(tj::BasisSymmetricNp<uint64_t>(
          nsites, nup, ndn, irrep.group(), characters));
    }
  } else {
    XDIAG_THROW(fmt::format("Unknown backend: \"{}\"", backend));
  }
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
tJ::iterator_t tJ::begin() const { return iterator_t(*this, true); }
tJ::iterator_t tJ::end() const { return iterator_t(*this, false); }
int64_t tJ::index(ProductState const &pstate) const try {
  return std::visit(
      [&](auto &&basis) {
        using basis_t = typename std::decay<decltype(basis)>::type;
        using bit_t = typename basis_t::bit_t;
        auto [ups, dns] = to_bits_tj<bit_t>(pstate);
        return basis.index(ups, dns);
      },
      *basis_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
bool tJ::operator==(tJ const &rhs) const {
  return (nsites_ == rhs.nsites_) && (nup_ == rhs.nup_) && (ndn_ == rhs.ndn_) &&
         (irrep_ == rhs.irrep_);
}
bool tJ::operator!=(tJ const &rhs) const { return !operator==(rhs); }
int64_t tJ::dim() const { return size_; }
int64_t tJ::size() const { return size_; }

int64_t tJ::nsites() const { return nsites_; }
std::string tJ::backend() const { return backend_; }
std::optional<int64_t> tJ::nup() const { return nup_; }
std::optional<int64_t> tJ::ndn() const { return ndn_; }
std::optional<Representation> const &tJ::irrep() const { return irrep_; }

bool tJ::isreal() const { return irrep_ ? irrep_->isreal() : true; }
tJ::basis_t const &tJ::basis() const { return *basis_; }

int64_t index(tJ const &block, ProductState const &pstate) {
  return block.index(pstate);
}
int64_t nsites(tJ const &block) { return block.nsites(); }
int64_t dim(tJ const &block) { return block.dim(); }
int64_t size(tJ const &block) { return block.size(); }
bool isreal(tJ const &block) { return block.isreal(); }

std::ostream &operator<<(std::ostream &out, tJ const &block) {
  out << "tJ:\n";
  out << "  nsites   : " << block.nsites() << "\n";
  if (block.nup()) {
    out << "  nup      : " << *block.nup() << "\n";
  } else {
    out << "  nup      : not conserved\n";
  }

  if (block.ndn()) {
    out << "  ndn      : " << *block.ndn() << "\n";
  } else {
    out << "  ndn      : not conserved\n";
  }
  if (block.irrep()) {
    out << "  irrep    : defined with ID " << std::hex
        << random::hash(*block.irrep()) << std::dec << "\n";
  }
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << block.size();
  out << "  dimension: " << ss.str() << "\n";
  out << "  ID       : " << std::hex << random::hash(block) << std::dec << "\n";
  return out;
}
std::string to_string(tJ const &block) { return to_string_generic(block); }

tJIterator::tJIterator(tJ const &block, bool begin)
    : nsites_(block.nsites()), pstate_(nsites_),
      it_(std::visit(
          [&](auto const &basis) {
            basis::BasistJIterator it = begin ? basis.begin() : basis.end();
            return it;
          },
          block.basis())) {}

tJIterator &tJIterator::operator++() {
  std::visit([](auto &&it) { ++it; }, it_);
  return *this;
}

ProductState const &tJIterator::operator*() const {
  std::visit(
      [&](auto &&it) {
        auto [ups, dns] = *it;
        to_product_state_tj(ups, dns, pstate_);
      },
      it_);
  return pstate_;
}

bool tJIterator::operator!=(tJIterator const &rhs) const {
  return it_ != rhs.it_;
}

} // namespace xdiag
