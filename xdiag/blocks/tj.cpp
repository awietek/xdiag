#include "tj.hpp"

#include <xdiag/random/hash.hpp>

namespace xdiag {

using namespace basis;

tJ::tJ(int64_t n_sites, int64_t nup, int64_t ndn) try
    : n_sites_(n_sites), n_up_(nup), n_dn_(ndn), irrep_(std::nullopt) {

  if (n_sites < 0) {
    XDIAG_THROW("Invalid argument: n_sites < 0");
  } else if ((nup < 0) || (ndn < 0)) {
    XDIAG_THROW("Invalid argument: (nup < 0) or (ndn < 0)");
  } else if ((nup + ndn) > n_sites) {
    XDIAG_THROW("Invalid argument: nup + ndn > n_sites");
  } else if (n_sites < 32) {
    basis_ =
        std::make_shared<basis_t>(tj::BasisNp<uint32_t>(n_sites, nup, ndn));
  } else if (n_sites < 64) {
    basis_ =
        std::make_shared<basis_t>(tj::BasisNp<uint64_t>(n_sites, nup, ndn));
  } else {
    XDIAG_THROW("blocks with more than 64 sites currently not implemented");
  }
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

tJ::tJ(int64_t n_sites, int64_t nup, int64_t ndn,
       Representation const &irrep) try
    : n_sites_(n_sites), n_up_(nup), n_dn_(ndn), irrep_(irrep) {

  if (n_sites < 0) {
    XDIAG_THROW("Invalid argument: n_sites < 0");
  } else if ((nup < 0) || (ndn < 0)) {
    XDIAG_THROW("Invalid argument: (nup < 0) or (ndn < 0)");
  } else if ((nup + ndn) > n_sites) {
    XDIAG_THROW("Invalid argument: nup + ndn > n_sites");
  } else if (n_sites != irrep.group().n_sites()) {
    XDIAG_THROW("n_sites does not match the n_sites in PermutationGroup");
  } else if (n_sites < 32) {
    if (irrep.isreal()) {
      auto characters = irrep.characters().as<arma::vec>();
      basis_ = std::make_shared<basis_t>(tj::BasisSymmetricNp<uint32_t>(
          n_sites, nup, ndn, irrep.group(), characters));
    } else {
      auto characters = irrep.characters().as<arma::cx_vec>();
      basis_ = std::make_shared<basis_t>(tj::BasisSymmetricNp<uint32_t>(
          n_sites, nup, ndn, irrep.group(), characters));
    }

  } else if (n_sites < 64) {
    if (irrep.isreal()) {
      auto characters = irrep.characters().as<arma::vec>();
      basis_ = std::make_shared<basis_t>(tj::BasisSymmetricNp<uint64_t>(
          n_sites, nup, ndn, irrep.group(), characters));
    } else {
      auto characters = irrep.characters().as<arma::cx_vec>();
      basis_ = std::make_shared<basis_t>(tj::BasisSymmetricNp<uint64_t>(
          n_sites, nup, ndn, irrep.group(), characters));
    }
  } else {
    XDIAG_THROW("blocks with more than 64 sites currently not implemented");
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
  return (n_sites_ == rhs.n_sites_) && (n_up_ == rhs.n_up_) &&
         (n_dn_ == rhs.n_dn_) && (irrep_ == rhs.irrep_);
}
bool tJ::operator!=(tJ const &rhs) const { return !operator==(rhs); }
int64_t tJ::dim() const { return size_; }
int64_t tJ::size() const { return size_; }

int64_t tJ::n_sites() const { return n_sites_; }
std::optional<int64_t> tJ::n_up() const { return n_up_; }
std::optional<int64_t> tJ::n_dn() const { return n_dn_; }
std::optional<Representation> const &tJ::irrep() const { return irrep_; }

bool tJ::isreal() const { return irrep_ ? irrep_->isreal() : true; }
tJ::basis_t const &tJ::basis() const { return *basis_; }

bool isreal(tJ const &block) { return block.isreal(); }

std::ostream &operator<<(std::ostream &out, tJ const &block) {
  out << "tJ:\n";
  out << "  n_sites  : " << block.n_sites() << "\n";
  if (block.n_up()) {
    out << "  n_up     : " << *block.n_up() << "\n";
  } else {
    out << "  n_up     : not conserved\n";
  }

  if (block.n_dn()) {
    out << "  n_dn     : " << *block.n_dn() << "\n";
  } else {
    out << "  n_dn     : not conserved\n";
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
    : n_sites_(block.n_sites()), pstate_(n_sites_),
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
