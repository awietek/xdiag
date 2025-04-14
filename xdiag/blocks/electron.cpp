// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "electron.hpp"
#include <xdiag/random/hash.hpp>

namespace xdiag {

using namespace basis;

Electron::Electron(int64_t nsites, std::string backend) try
    : nsites_(nsites), backend_(backend), nup_(std::nullopt),
      ndn_(std::nullopt), irrep_(std::nullopt) {
  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  }

  // Choose basis implementation
  if (backend == "auto") {
    if (nsites < 32) {
      basis_ = std::make_shared<basis_t>(electron::BasisNoNp<uint32_t>(nsites));
    } else if (nsites < 64) {
      basis_ = std::make_shared<basis_t>(electron::BasisNoNp<uint64_t>(nsites));
    } else {
      XDIAG_THROW(
          "Spinhalf blocks with more than 64 sites currently not implemented");
    }
  } else if (backend == "32bit") {
    basis_ = std::make_shared<basis_t>(electron::BasisNoNp<uint32_t>(nsites));
  } else if (backend == "64bit") {
    basis_ = std::make_shared<basis_t>(electron::BasisNoNp<uint64_t>(nsites));
  } else {
    XDIAG_THROW(fmt::format("Unknown backend: \"{}\"", backend));
  }
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Electron::Electron(int64_t nsites, int64_t nup, int64_t ndn,
                   std::string backend) try
    : nsites_(nsites), backend_(backend), nup_(nup), ndn_(ndn),
      irrep_(std::nullopt) {

  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  } else if ((nup < 0) || (nup > nsites)) {
    XDIAG_THROW("Invalid argument: (nup < 0) or (nup > nsites)");
  } else if ((ndn < 0) || (ndn > nsites)) {
    XDIAG_THROW("Invalid argument: (ndn < 0) or (ndn > nsites)");
  }

  // Choose basis implementation
  if (backend == "auto") {
    if (nsites < 32) {
      basis_ = std::make_shared<basis_t>(
          electron::BasisNp<uint32_t>(nsites, nup, ndn));
    } else if (nsites < 64) {
      basis_ = std::make_shared<basis_t>(
          electron::BasisNp<uint64_t>(nsites, nup, ndn));
    } else {
      XDIAG_THROW("Blocks with more than 64 sites currently not implemented");
    }
  } else if (backend == "32bit") {
    basis_ = std::make_shared<basis_t>(
        electron::BasisNp<uint32_t>(nsites, nup, ndn));
  } else if (backend == "64bit") {
    basis_ = std::make_shared<basis_t>(
        electron::BasisNp<uint64_t>(nsites, nup, ndn));
  } else {
    XDIAG_THROW(fmt::format("Unknown backend: \"{}\"", backend));
  }
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Electron::Electron(int64_t nsites, Representation const &irrep,
                   std::string backend) try
    : nsites_(nsites), backend_(backend), nup_(std::nullopt),
      ndn_(std::nullopt), irrep_(irrep) {
  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  } else if (nsites != irrep.group().nsites()) {
    XDIAG_THROW("nsites does not match the nsites in PermutationGroup");
  }

  // Choose basis implementation
  if (backend == "auto") {
    if (nsites < 32) {
      basis_ = std::make_shared<basis_t>(
          electron::BasisSymmetricNoNp<uint32_t>(nsites, irrep));
    } else if (nsites < 64) {
      basis_ = std::make_shared<basis_t>(
          electron::BasisSymmetricNoNp<uint64_t>(nsites, irrep));
    } else {
      XDIAG_THROW("Blocks with more than 64 sites currently not implemented");
    }
  } else if (backend == "32bit") {
    basis_ = std::make_shared<basis_t>(
        electron::BasisSymmetricNoNp<uint32_t>(nsites, irrep));
  } else if (backend == "64bit") {
    basis_ = std::make_shared<basis_t>(
        electron::BasisSymmetricNoNp<uint64_t>(nsites, irrep));
  } else {
    XDIAG_THROW(fmt::format("Unknown backend: \"{}\"", backend));
  }
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Electron::Electron(int64_t nsites, int64_t nup, int64_t ndn,
                   Representation const &irrep, std::string backend) try
    : nsites_(nsites), backend_(backend), nup_(nup), ndn_(ndn), irrep_(irrep) {
  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  } else if ((nup < 0) || (nup > nsites)) {
    XDIAG_THROW("Invalid argument: (nup < 0) or (nup > nsites)");
  } else if ((ndn < 0) || (ndn > nsites)) {
    XDIAG_THROW("Invalid argument: (ndn < 0) or (ndn > nsites)");
  } else if (nsites != irrep.group().nsites()) {
    XDIAG_THROW("nsites does not match the nsites in PermutationGroup");
  }

  // Choose basis implementation
  if (backend == "auto") {
    if (nsites < 32) {
      basis_ = std::make_shared<basis_t>(
          electron::BasisSymmetricNp<uint32_t>(nsites, nup, ndn, irrep));
    } else if (nsites < 64) {
      basis_ = std::make_shared<basis_t>(
          electron::BasisSymmetricNp<uint64_t>(nsites, nup, ndn, irrep));
    } else {
      XDIAG_THROW("Blocks with more than 64 sites currently not implemented");
    }
  } else if (backend == "32bit") {
    basis_ = std::make_shared<basis_t>(
        electron::BasisSymmetricNp<uint32_t>(nsites, nup, ndn, irrep));
  } else if (backend == "64bit") {
    basis_ = std::make_shared<basis_t>(
        electron::BasisSymmetricNp<uint64_t>(nsites, nup, ndn, irrep));
  } else {
    XDIAG_THROW(fmt::format("Unknown backend: \"{}\"", backend));
  }
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Electron::iterator_t Electron::begin() const { return iterator_t(*this, true); }
Electron::iterator_t Electron::end() const { return iterator_t(*this, false); }
int64_t Electron::index(ProductState const &pstate) const try {
  return std::visit(
      [&](auto &&basis) {
        using basis_t = typename std::decay<decltype(basis)>::type;
        using bit_t = typename basis_t::bit_t;
        auto [ups, dns] = to_bits_electron<bit_t>(pstate);
        return basis.index(ups, dns);
      },
      *basis_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool Electron::operator==(Electron const &rhs) const {
  return (nsites_ == rhs.nsites_) && (nup_ == rhs.nup_) && (ndn_ == rhs.ndn_) &&
         (irrep_ == rhs.irrep_);
}
bool Electron::operator!=(Electron const &rhs) const {
  return !operator==(rhs);
}
int64_t Electron::dim() const { return size_; }
int64_t Electron::size() const { return size_; }

int64_t Electron::nsites() const { return nsites_; }
std::string Electron::backend() const { return backend_; }
std::optional<int64_t> Electron::nup() const { return nup_; }
std::optional<int64_t> Electron::ndn() const { return ndn_; }
std::optional<Representation> const &Electron::irrep() const { return irrep_; }
bool Electron::isreal() const { return irrep_ ? irrep_->isreal() : true; }
Electron::basis_t const &Electron::basis() const { return *basis_; }

int64_t index(Electron const &block, ProductState const &pstate) {
  return block.index(pstate);
}
int64_t nsites(Electron const &block) { return block.nsites(); }
int64_t dim(Electron const &block) { return block.dim(); }
int64_t size(Electron const &block) { return block.size(); }
bool isreal(Electron const &block) { return block.isreal(); }
std::ostream &operator<<(std::ostream &out, Electron const &block) {
  out << "Electron:\n";
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
std::string to_string(Electron const &block) {
  return to_string_generic(block);
}

ElectronIterator::ElectronIterator(Electron const &block, bool begin)
    : nsites_(block.nsites()), pstate_(nsites_),
      it_(std::visit(
          [&](auto const &basis) {
            basis::BasisElectronIterator it =
                begin ? basis.begin() : basis.end();
            return it;
          },
          block.basis())) {}

ElectronIterator &ElectronIterator::operator++() {
  std::visit([](auto &&it) { ++it; }, it_);
  return *this;
}

ProductState const &ElectronIterator::operator*() const {
  std::visit(
      [&](auto &&it) {
        auto [ups, dns] = *it;
        to_product_state_electron(ups, dns, pstate_);
      },
      it_);
  return pstate_;
}

bool ElectronIterator::operator!=(ElectronIterator const &rhs) const {
  return it_ != rhs.it_;
}

} // namespace xdiag
