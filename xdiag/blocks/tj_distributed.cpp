#include "tj_distributed.hpp"

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/random/hash.hpp>

namespace xdiag {

tJDistributed::tJDistributed(int64_t nsites, int64_t nup, int64_t ndn,
                             std::string backend) try
    : nsites_(nsites), backend_(backend), nup_(nup), ndn_(ndn) {
  using namespace basis::tj_distributed;
  using combinatorics::binomial;

  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("nsites < 0");
  } else if ((nup < 0) || (ndn < 0)) {
    XDIAG_THROW("nup < 0 or ndn < 0");
  } else if ((nup + ndn) > nsites) {
    XDIAG_THROW("nup + ndn > nsites");
  }
  // Choose basis implementation
  if (backend == "auto") {
    if (nsites < 32) {
      basis_ =
          std::make_shared<basis_t>(BasisNp<uint32_t>(nsites, nup, ndn));
    } else if (nsites < 64) {
      basis_ =
          std::make_shared<basis_t>(BasisNp<uint64_t>(nsites, nup, ndn));
    } else {
      XDIAG_THROW("Blocks with more than 64 sites currently not implemented");
    }
  } else if (backend == "32bit") {
    basis_ = std::make_shared<basis_t>(BasisNp<uint32_t>(nsites, nup, ndn));
  } else if (backend == "64bit") {
    basis_ = std::make_shared<basis_t>(BasisNp<uint64_t>(nsites, nup, ndn));
  } else {
    XDIAG_THROW(fmt::format("Unknown backend: \"{}\"", backend));
  }
  dim_ = basis::dim(*basis_);
  assert(dim_ == binomial(nsites, nup) * binomial(nsites - nup, ndn));
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t tJDistributed::nsites() const { return nsites_; }
int64_t tJDistributed::nup() const { return nup_; }
int64_t tJDistributed::ndn() const { return ndn_; }

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
  return (nsites_ == rhs.nsites_) && (nup_ == rhs.nup_) &&
         (ndn_ == rhs.ndn_);
}
bool tJDistributed::operator!=(tJDistributed const &rhs) const {
  return !operator==(rhs);
}

tJDistributed::basis_t const &tJDistributed::basis() const { return *basis_; }

bool isreal(tJDistributed const &block) { return block.isreal(); }
std::ostream &operator<<(std::ostream &out, tJDistributed const &block) {
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  out << "tJDistributed:\n";
  out << "  nsites  : " << block.nsites() << "\n";
  if ((block.nup() != undefined) && (block.ndn() != undefined)) {
    out << "  nup     : " << block.nup() << "\n";
    out << "  ndn     : " << block.ndn() << "\n";
  } else {
    out << "  nup     : not conserved\n";
    out << "  ndn     : not conserved\n";
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
    : nsites_(block.nsites()), pstate_(nsites_),
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
