#include "electron.hpp"
#include <xdiag/random/hash.hpp>

namespace xdiag {

using namespace basis;

Electron::Electron(int64_t n_sites) try
    : n_sites_(n_sites), n_up_(undefined), n_dn_(undefined),
      permutation_group_(), irrep_() {

  if (n_sites < 0) {
    XDIAG_THROW("Invalid argument: n_sites < 0");
  } else if (n_sites < 32) {
    basis_ = std::make_shared<basis_t>(electron::BasisNoNp<uint32_t>(n_sites));
  } else if (n_sites < 64) {
    basis_ = std::make_shared<basis_t>(electron::BasisNoNp<uint64_t>(n_sites));
  } else {
    XDIAG_THROW(
        "Spinhalf blocks with more than 64 sites currently not implemented");
  }
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Electron::Electron(int64_t n_sites, int64_t nup, int64_t ndn) try
    : n_sites_(n_sites), n_up_(nup), n_dn_(ndn), permutation_group_(),
      irrep_() {

  if (n_sites < 0) {
    XDIAG_THROW("Invalid argument: n_sites < 0");
  } else if ((nup < 0) || (nup > n_sites)) {
    XDIAG_THROW("Invalid argument: (nup < 0) or (nup > n_sites)");
  } else if ((ndn < 0) || (ndn > n_sites)) {
    XDIAG_THROW("Invalid argument: (ndn < 0) or (ndn > n_sites)");
  } else if (n_sites < 32) {
    basis_ = std::make_shared<basis_t>(
        electron::BasisNp<uint32_t>(n_sites, nup, ndn));
  } else if (n_sites < 64) {
    basis_ = std::make_shared<basis_t>(
        electron::BasisNp<uint64_t>(n_sites, nup, ndn));
  } else {
    XDIAG_THROW("blocks with more than 64 sites currently not implemented");
  }
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Electron::Electron(int64_t n_sites, PermutationGroup group,
                   Representation irrep) try
    : n_sites_(n_sites), n_up_(undefined), n_dn_(undefined),
      permutation_group_(allowed_subgroup(group, irrep)), irrep_(irrep) {

  if (n_sites < 0) {
    XDIAG_THROW("Invalid argument: n_sites < 0");
  } else if (n_sites != group.n_sites()) {
    XDIAG_THROW("n_sites does not match the n_sites in PermutationGroup");
  } else if (permutation_group_.size() != irrep.size()) {
    XDIAG_THROW("PermutationGroup and Representation do not have "
                "same number of elements");
  } else if (n_sites < 32) {
    basis_ = std::make_shared<basis_t>(
        electron::BasisSymmetricNoNp<uint32_t>(n_sites, group, irrep));
  } else if (n_sites < 64) {
    basis_ = std::make_shared<basis_t>(
        electron::BasisSymmetricNoNp<uint64_t>(n_sites, group, irrep));
  } else {
    XDIAG_THROW("blocks with more than 64 sites currently not implemented");
  }
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Electron::Electron(int64_t n_sites, int64_t nup, int64_t ndn,
                   PermutationGroup group, Representation irrep) try
    : n_sites_(n_sites), n_up_(nup), n_dn_(ndn),
      permutation_group_(allowed_subgroup(group, irrep)), irrep_(irrep) {

  if (n_sites < 0) {
    XDIAG_THROW("Invalid argument: n_sites < 0");
  } else if ((nup < 0) || (nup > n_sites)) {
    XDIAG_THROW("Invalid argument: (nup < 0) or (nup > n_sites)");
  } else if ((ndn < 0) || (ndn > n_sites)) {
    XDIAG_THROW("Invalid argument: (ndn < 0) or (ndn > n_sites)");
  } else if (n_sites != group.n_sites()) {
    XDIAG_THROW("n_sites does not match the n_sites in PermutationGroup");
  } else if (permutation_group_.size() != irrep.size()) {
    XDIAG_THROW("PermutationGroup and Representation do not have "
                "same number of elements");
  } else if (n_sites < 32) {
    basis_ = std::make_shared<basis_t>(
        electron::BasisSymmetricNp<uint32_t>(n_sites, nup, ndn, group, irrep));
  } else if (n_sites < 64) {
    basis_ = std::make_shared<basis_t>(
        electron::BasisSymmetricNp<uint64_t>(n_sites, nup, ndn, group, irrep));
  } else {
    XDIAG_THROW("blocks with more than 64 sites currently not implemented");
  }
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t Electron::n_sites() const { return n_sites_; }
int64_t Electron::n_up() const { return n_up_; }
int64_t Electron::n_dn() const { return n_dn_; }

PermutationGroup const &Electron::permutation_group() const {
  return permutation_group_;
}
Representation const &Electron::irrep() const { return irrep_; }

int64_t Electron::dim() const { return size_; }
int64_t Electron::size() const { return size_; }
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
bool Electron::isreal(double precision) const {
  return irrep_.isreal(precision);
}

bool Electron::operator==(Electron const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_up_ == rhs.n_up_) &&
         (n_dn_ == rhs.n_dn_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}
bool Electron::operator!=(Electron const &rhs) const {
  return !operator==(rhs);
}

Electron::basis_t const &Electron::basis() const { return *basis_; }

std::ostream &operator<<(std::ostream &out, Electron const &block) {
  out << "Electron:\n";
  out << "  n_sites  : " << block.n_sites() << "\n";
  if ((block.n_up() != undefined) && (block.n_dn() != undefined)) {
    out << "  n_up     : " << block.n_up() << "\n";
    out << "  n_dn     : " << block.n_dn() << "\n";
  } else {
    out << "  n_up     : not conserved\n";
    out << "  n_dn     : not conserved\n";
  }
  if (block.permutation_group()) {
    out << "  group    : defined with ID " << std::hex
        << random::hash(block.permutation_group()) << std::dec << "\n";
    out << "  irrep    : defined with ID " << std::hex
        << random::hash(block.irrep()) << std::dec << "\n";
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
    : n_sites_(block.n_sites()), pstate_(n_sites_),
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
