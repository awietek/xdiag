#include "tj_distributed.hpp"

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/random/hash.hpp>

namespace xdiag {

tJDistributed::tJDistributed(int64_t n_sites, int64_t n_up, int64_t n_dn) try
    : n_sites_(n_sites), n_up_(n_up), n_dn_(n_dn) {
  using namespace basis::tj_distributed;
  using combinatorics::binomial;

  if (n_sites < 0) {
    XDIAG_THROW("n_sites < 0");
  } else if ((n_up < 0) || (n_dn < 0)) {
    XDIAG_THROW("n_up < 0 or n_dn < 0");
  } else if ((n_up + n_dn) > n_sites) {
    XDIAG_THROW("n_up + n_dn > n_sites");
  }

  if (n_sites < 32) {
    basis_ = std::make_shared<basis_t>(BasisNp<uint32_t>(n_sites, n_up, n_dn));
  } else if (n_sites < 64) {
    basis_ = std::make_shared<basis_t>(BasisNp<uint64_t>(n_sites, n_up, n_dn));
  } else {
    XDIAG_THROW("blocks with more than 64 sites currently not implemented");
  }
  dim_ = basis::dim(*basis_);
  assert(dim_ == binomial(n_sites, n_up) * binomial(n_sites - n_up, n_dn));
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t tJDistributed::n_sites() const { return n_sites_; }
int64_t tJDistributed::n_up() const { return n_up_; }
int64_t tJDistributed::n_dn() const { return n_dn_; }

int64_t tJDistributed::dim() const { return dim_; }
int64_t tJDistributed::size() const { return size_; }
int64_t tJDistributed::size_max() const { return basis::size_max(*basis_); }
int64_t tJDistributed::size_min() const { return basis::size_min(*basis_); }
tJDistributed::iterator_t tJDistributed::begin() const {
  return iterator_t(*this, true);
}
tJDistributed::iterator_t tJDistributed::end() const {
  return iterator_t(*this, false);
}
int64_t tJDistributed::index(ProductState const &pstate) const try {
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
bool tJDistributed::isreal(double precision) const {
  return true; // would only be nontrivial with space group irreps
}

bool tJDistributed::operator==(tJDistributed const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_up_ == rhs.n_up_) &&
         (n_dn_ == rhs.n_dn_);
}
bool tJDistributed::operator!=(tJDistributed const &rhs) const {
  return !operator==(rhs);
}

tJDistributed::basis_t const &tJDistributed::basis() const { return *basis_; }

std::ostream &operator<<(std::ostream &out, tJDistributed const &block) {
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  out << "tJDistributed:\n";
  out << "  n_sites  : " << block.n_sites() << "\n";
  if ((block.n_up() != undefined) && (block.n_dn() != undefined)) {
    out << "  n_up     : " << block.n_up() << "\n";
    out << "  n_dn     : " << block.n_dn() << "\n";
  } else {
    out << "  n_up     : not conserved\n";
    out << "  n_dn     : not conserved\n";
  }

  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << block.dim();

  std::stringstream ssmax;
  ssmax.imbue(std::locale("en_US.UTF-8"));
  ssmax << block.size_max();

  std::stringstream ssmin;
  ssmin.imbue(std::locale("en_US.UTF-8"));
  ssmin << block.size_min();

  std::stringstream ssavg;
  ssavg.imbue(std::locale("en_US.UTF-8"));
  ssavg << block.dim() / mpi_size;

  out << "  dimension       : " << ss.str() << "\n";
  out << "  size (max local): " << ssmax.str() << "\n";
  out << "  size (min local): " << ssmin.str() << "\n";
  out << "  size (avg local): " << ssavg.str() << "\n";
  out << "  ID              : " << std::hex << random::hash(block) << std::dec
      << "\n";
  return out;
}
std::string to_string(tJDistributed const &block) {
  return to_string_generic(block);
}

tJDistributedIterator::tJDistributedIterator(tJDistributed const &block,
                                             bool begin)
    : n_sites_(block.n_sites()), pstate_(n_sites_),
      it_(std::visit(
          [&](auto const &basis) {
            basis::BasistJDistributedIterator it =
                begin ? basis.begin() : basis.end();
            return it;
          },
          block.basis())) {}

tJDistributedIterator &tJDistributedIterator::operator++() {
  std::visit([](auto &&it) { ++it; }, it_);
  return *this;
}

ProductState const &tJDistributedIterator::operator*() const {
  std::visit(
      [&](auto &&it) {
        auto [ups, dns] = *it;
        to_product_state_tj(ups, dns, pstate_);
      },
      it_);
  return pstate_;
}

bool tJDistributedIterator::operator!=(tJDistributedIterator const &rhs) const {
  return it_ != rhs.it_;
}

} // namespace xdiag
