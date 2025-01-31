#include "spinhalf_distributed.hpp"
#include <xdiag/random/hash.hpp>

namespace xdiag {

SpinhalfDistributed::SpinhalfDistributed(int64_t nsites, int64_t nup,
                                         std::string backend) try
    : nsites_(nsites), backend_(backend), nup_(nup) {
  using namespace basis::spinhalf_distributed;
  using combinatorics::binomial;

  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("nsites < 0");
  } else if (nup < 0) {
    XDIAG_THROW("nup < 0");
  } else if (nup > nsites) {
    XDIAG_THROW("nup > nsites");
  }
  // Choose basis implementation
  if (backend == "auto") {
    if (nsites < 32) {
      basis_ = std::make_shared<basis_t>(BasisSz<uint32_t>(nsites, nup));
    } else if (nsites < 64) {
      basis_ = std::make_shared<basis_t>(BasisSz<uint64_t>(nsites, nup));
    } else {
      XDIAG_THROW("Blocks with more than 64 sites currently not implemented");
    }
  } else if (backend == "32bit") {
    basis_ = std::make_shared<basis_t>(BasisSz<uint32_t>(nsites, nup));
  } else if (backend == "64bit") {
    basis_ = std::make_shared<basis_t>(BasisSz<uint64_t>(nsites, nup));
  } else {
    XDIAG_THROW(fmt::format("Unknown backend: \"{}\"", backend));
  }

  dim_ = basis::dim(*basis_);
  assert(dim_ == binomial(nsites, nup));
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t SpinhalfDistributed::nsites() const { return nsites_; }
std::string SpinhalfDistributed::backend() const { return backend_; }
std::optional<int64_t> SpinhalfDistributed::nup() const { return nup_; }

int64_t SpinhalfDistributed::dim() const { return dim_; }
int64_t SpinhalfDistributed::size() const { return size_; }
int64_t SpinhalfDistributed::size_max() const {
  return basis::size_max(*basis_);
}
int64_t SpinhalfDistributed::size_min() const {
  return basis::size_min(*basis_);
}
SpinhalfDistributed::iterator_t SpinhalfDistributed::begin() const {
  return iterator_t(*this, true);
}
SpinhalfDistributed::iterator_t SpinhalfDistributed::end() const {
  return iterator_t(*this, false);
}
int64_t SpinhalfDistributed::index(ProductState const &pstate) const try {
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
bool SpinhalfDistributed::isreal(double precision) const {
  return true; // would only be nontrivial with space group irreps
}

bool SpinhalfDistributed::operator==(SpinhalfDistributed const &rhs) const {
  return (nsites_ == rhs.nsites_) && (nup_ == rhs.nup_);
}
bool SpinhalfDistributed::operator!=(SpinhalfDistributed const &rhs) const {
  return !operator==(rhs);
}

SpinhalfDistributed::basis_t const &SpinhalfDistributed::basis() const {
  return *basis_;
}

int64_t index(SpinhalfDistributed const &block, ProductState const &pstate) {
  return block.index(pstate);
}
int64_t nsites(SpinhalfDistributed const &block) { return block.nsites(); }
int64_t dim(SpinhalfDistributed const &block) { return block.dim(); }
int64_t size(SpinhalfDistributed const &block) { return block.size(); }
bool isreal(SpinhalfDistributed const &block) { return block.isreal(); }
std::ostream &operator<<(std::ostream &out, SpinhalfDistributed const &block) {
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  out << "SpinhalfDistributed:\n";
  out << "  nsites   : " << block.nsites() << "\n";
  if (block.nup()) {
    out << "  nup      : " << *block.nup() << "\n";
  } else {
    out << "  nup      : not conserved\n";
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
std::string to_string(SpinhalfDistributed const &block) {
  return to_string_generic(block);
}

SpinhalfDistributedIterator::SpinhalfDistributedIterator(
    SpinhalfDistributed const &block, bool begin)
    : nsites_(block.nsites()), pstate_(nsites_),
      it_(std::visit(
          [&](auto const &basis) {
            basis::BasisSpinhalfDistributedIterator it =
                begin ? basis.begin() : basis.end();
            return it;
          },
          block.basis())) {}

SpinhalfDistributedIterator &SpinhalfDistributedIterator::operator++() {
  std::visit([](auto &&it) { ++it; }, it_);
  return *this;
}

ProductState const &SpinhalfDistributedIterator::operator*() const {
  std::visit(
      [&](auto &&it) {
        auto spins = *it;
        to_product_state_spinhalf(spins, pstate_);
      },
      it_);
  return pstate_;
}

bool SpinhalfDistributedIterator::operator!=(
    SpinhalfDistributedIterator const &rhs) const {
  return it_ != rhs.it_;
}

} // namespace xdiag
