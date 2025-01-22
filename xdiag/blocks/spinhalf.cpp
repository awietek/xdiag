#include "spinhalf.hpp"

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/random/hash.hpp>

namespace xdiag {

using namespace basis;

Spinhalf::Spinhalf(int64_t nsites, std::string backend) try
    : nsites_(nsites), backend_(backend), nup_(std::nullopt),
      irrep_(std::nullopt), size_((int64_t)1 << nsites) {
  check_dimension_works_with_blas_int_size(size_);

  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  }

  // Choose basis implementation
  if (backend == "auto") {
    if (nsites < 32) {
      basis_ = std::make_shared<basis_t>(spinhalf::BasisNoSz<uint32_t>(nsites));
    } else if (nsites < 64) {
      basis_ = std::make_shared<basis_t>(spinhalf::BasisNoSz<uint64_t>(nsites));
    } else {
      XDIAG_THROW(
          "Spinhalf blocks with more than 64 sites currently not implemented");
    }
  } else if (backend == "32bit") {
    basis_ = std::make_shared<basis_t>(spinhalf::BasisNoSz<uint32_t>(nsites));
  } else if (backend == "64bit") {
    basis_ = std::make_shared<basis_t>(spinhalf::BasisNoSz<uint64_t>(nsites));
  } else {
    XDIAG_THROW(fmt::format("Unknown backend: \"{}\"", backend));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Spinhalf::Spinhalf(int64_t nsites, int64_t nup, std::string backend) try
    : nsites_(nsites), backend_(backend), nup_(nup), irrep_(std::nullopt),
      size_(combinatorics::binomial(nsites, nup)) {
  check_dimension_works_with_blas_int_size(size_);

  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  } else if (nup < 0) {
    XDIAG_THROW("Invalid argument: nup < 0");
  } else if (nup > nsites) {
    XDIAG_THROW("Invalid argument: nup > nsites");
  }

  // Choose basis implementation
  if (backend == "auto") {
    if (nsites < 32) {
      basis_ =
          std::make_shared<basis_t>(spinhalf::BasisSz<uint32_t>(nsites, nup));
    } else if (nsites < 64) {
      basis_ =
          std::make_shared<basis_t>(spinhalf::BasisSz<uint64_t>(nsites, nup));
    } else {
      XDIAG_THROW(
          "Spinhalf blocks with more than 64 sites currently not implemented");
    }
  } else if (backend == "32bit") {
    basis_ =
        std::make_shared<basis_t>(spinhalf::BasisSz<uint32_t>(nsites, nup));
  } else if (backend == "64bit") {
    basis_ =
        std::make_shared<basis_t>(spinhalf::BasisSz<uint64_t>(nsites, nup));
  } else {
    XDIAG_THROW(fmt::format("Unknown backend: \"{}\"", backend));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Spinhalf::Spinhalf(int64_t nsites, Representation const &irrep,
                   std::string backend) try
    : nsites_(nsites), backend_(backend), nup_(std::nullopt), irrep_(irrep) {

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
          spinhalf::BasisSymmetricNoSz<uint32_t>(irrep));
    } else if (nsites < 64) {
      basis_ = std::make_shared<basis_t>(
          spinhalf::BasisSymmetricNoSz<uint64_t>(irrep));
    } else {
      XDIAG_THROW(
          "Spinhalf blocks with more than 64 sites currently not implemented");
    }
  } else if (backend == "32bit") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSymmetricNoSz<uint32_t>(irrep));
  } else if (backend == "64bit") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSymmetricNoSz<uint64_t>(irrep));
  } else if (backend == "1sublattice") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSublattice<uint64_t, 1>(irrep));
  } else if (backend == "2sublattice") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSublattice<uint64_t, 2>(irrep));
  } else if (backend == "3sublattice") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSublattice<uint64_t, 3>(irrep));
  } else if (backend == "4sublattice") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSublattice<uint64_t, 4>(irrep));
  } else if (backend == "5sublattice") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSublattice<uint64_t, 5>(irrep));
  } else {
    XDIAG_THROW(fmt::format("Unknown backend: \"{}\"", backend));
  }

  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Spinhalf::Spinhalf(int64_t nsites, int64_t nup, Representation const &irrep,
                   std::string backend) try
    : nsites_(nsites), backend_(backend), nup_(nup), irrep_(irrep) {

  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("Invalid argument: nsites < 0");
  } else if (nup < 0) {
    XDIAG_THROW("Invalid argument: nup < 0");
  } else if (nup > nsites) {
    XDIAG_THROW("Invalid argument: nup > nsites");
  } else if (nsites != irrep.group().nsites()) {
    XDIAG_THROW("nsites does not match the nsites in PermutationGroup");
  }

  // Choose basis implementation
  if (backend == "auto") {
    if (nsites < 32) {
      basis_ = std::make_shared<basis_t>(
          spinhalf::BasisSymmetricSz<uint32_t>(nup, irrep));
    } else if (nsites < 64) {
      basis_ = std::make_shared<basis_t>(
          spinhalf::BasisSymmetricSz<uint64_t>(nup, irrep));
    } else {
      XDIAG_THROW(
          "Spinhalf blocks with more than 64 sites currently not implemented");
    }
  } else if (backend == "32bit") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSymmetricSz<uint32_t>(nup, irrep));
  } else if (backend == "64bit") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSymmetricSz<uint64_t>(nup, irrep));
  } else if (backend == "1sublattice") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSublattice<uint64_t, 1>(nup, irrep));
  } else if (backend == "2sublattice") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSublattice<uint64_t, 2>(nup, irrep));
  } else if (backend == "3sublattice") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSublattice<uint64_t, 3>(nup, irrep));
  } else if (backend == "4sublattice") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSublattice<uint64_t, 4>(nup, irrep));
  } else if (backend == "5sublattice") {
    basis_ = std::make_shared<basis_t>(
        spinhalf::BasisSublattice<uint64_t, 5>(nup, irrep));
  } else {
    XDIAG_THROW(fmt::format("Unknown backend: \"{}\"", backend));
  }

  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Spinhalf::iterator_t Spinhalf::begin() const { return iterator_t(*this, true); }
Spinhalf::iterator_t Spinhalf::end() const { return iterator_t(*this, false); }
int64_t Spinhalf::index(ProductState const &pstate) const try {
  return std::visit(
      [&](auto &&basis) {
        using basis_t = typename std::decay<decltype(basis)>::type;
        using bit_t = typename basis_t::bit_t;
        bit_t spins = to_bits_spinhalf<bit_t>(pstate);
        return basis.index(spins);
      },
      *basis_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
int64_t Spinhalf::dim() const { return size_; }
int64_t Spinhalf::size() const { return size_; }

bool Spinhalf::operator==(Spinhalf const &rhs) const {
  return (nsites_ == rhs.nsites_) && (nup_ == rhs.nup_) &&
         (irrep_ == rhs.irrep_) && (*basis_ == *rhs.basis_);
}

bool Spinhalf::operator!=(Spinhalf const &rhs) const {
  return !operator==(rhs);
}

int64_t Spinhalf::nsites() const { return nsites_; }
std::string Spinhalf::backend() const { return backend_; }
std::optional<int64_t> Spinhalf::nup() const { return nup_; }
std::optional<Representation> const &Spinhalf::irrep() const { return irrep_; }

int64_t index(Spinhalf const &block, ProductState const &pstate) {
  return block.index(pstate);
}
int64_t nsites(Spinhalf const &block) { return block.nsites(); }

int64_t dim(Spinhalf const &block) { return block.dim(); }
int64_t size(Spinhalf const &block) { return block.size(); }
bool Spinhalf::isreal() const { return irrep_ ? irrep_->isreal() : true; }
Spinhalf::basis_t const &Spinhalf::basis() const { return *basis_; }

bool isreal(Spinhalf const &block) { return block.isreal(); }
std::ostream &operator<<(std::ostream &out, Spinhalf const &block) {
  out << "Spinhalf:\n";
  out << "  nsites  : " << block.nsites() << "\n";
  if (block.nup()) {
    out << "  nup     : " << *block.nup() << "\n";
  } else {
    out << "  nup     : not conserved\n";
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
std::string to_string(Spinhalf const &block) {
  return to_string_generic(block);
}

SpinhalfIterator::SpinhalfIterator(Spinhalf const &block, bool begin)
    : nsites_(block.nsites()), pstate_(nsites_),
      it_(std::visit(
          [&](auto const &basis) {
            basis::BasisSpinhalfIterator it =
                begin ? basis.begin() : basis.end();
            return it;
          },
          block.basis())) {}

SpinhalfIterator &SpinhalfIterator::operator++() {
  std::visit([](auto &&it) { ++it; }, it_);
  return *this;
}

ProductState const &SpinhalfIterator::operator*() const {
  std::visit(
      [&](auto &&it) {
        auto spins = *it;
        to_product_state_spinhalf(spins, pstate_);
      },
      it_);
  return pstate_;
}

bool SpinhalfIterator::operator!=(SpinhalfIterator const &rhs) const {
  return it_ != rhs.it_;
}

} // namespace xdiag
