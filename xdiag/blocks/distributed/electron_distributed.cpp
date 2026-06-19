// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "electron_distributed.hpp"

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/random/hash.hpp>

namespace xdiag {

ElectronDistributed::ElectronDistributed(int64_t nsites, int64_t nup,
                                         int64_t ndn, std::string backend) try
    : nsites_(nsites), backend_(backend), nup_(nup), ndn_(ndn) {
  using namespace combinatorics;
  using namespace basis::electron_distributed;

  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("nsites < 0");
  } else if ((nup < 0) || (ndn < 0)) {
    XDIAG_THROW("nup < 0 or ndn < 0");
  } else if (nup > nsites) {
    XDIAG_THROW("nup > nsites");
  } else if (ndn > nsites) {
    XDIAG_THROW("ndn > nsites");
  }

  // Choose basis implementation
  if (backend == "auto") {
    if (nsites < 32) {
      basis_ = std::make_shared<basis_t>(BasisNp<uint32_t>(nsites, nup, ndn));
    } else if (nsites < 64) {
      basis_ = std::make_shared<basis_t>(BasisNp<uint64_t>(nsites, nup, ndn));
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
  assert(dim_ == binomial(nsites, nup) * binomial(nsites, ndn));
  size_ = basis::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
}
XDIAG_CATCH

int64_t ElectronDistributed::nsites() const { return nsites_; }
std::string ElectronDistributed::backend() const { return backend_; }
std::optional<int64_t> ElectronDistributed::nup() const { return nup_; }
std::optional<int64_t> ElectronDistributed::ndn() const { return ndn_; }

int64_t ElectronDistributed::dim() const { return dim_; }
int64_t ElectronDistributed::size() const { return size_; }
int64_t ElectronDistributed::size_max() const {
  return basis::size_max(*basis_);
}
int64_t ElectronDistributed::size_min() const {
  return basis::size_min(*basis_);
}
ElectronDistributed::iterator_t ElectronDistributed::begin() const {
  return iterator_t(*this, true);
}
ElectronDistributed::iterator_t ElectronDistributed::end() const {
  return iterator_t(*this, false);
}
int64_t ElectronDistributed::index(ProductState const &pstate) const try {
  return std::visit(
      [&](auto &&basis) {
        using basis_t = typename std::decay<decltype(basis)>::type;
        using bit_t = typename basis_t::bit_t;
        auto [ups, dns] = to_bits_electron<bit_t>(pstate);
        return basis.index(ups, dns);
      },
      *basis_);
}
XDIAG_CATCH

bool ElectronDistributed::isreal() const {
  return true; // would only be nontrivial with space group irreps
}

bool ElectronDistributed::operator==(ElectronDistributed const &rhs) const {
  return (nsites_ == rhs.nsites_) && (nup_ == rhs.nup_) && (ndn_ == rhs.ndn_);
}
bool ElectronDistributed::operator!=(ElectronDistributed const &rhs) const {
  return !operator==(rhs);
}

ElectronDistributed::basis_t const &ElectronDistributed::basis() const {
  return *basis_;
}
int64_t index(ElectronDistributed const &block, ProductState const &pstate) {
  return block.index(pstate);
}
int64_t nsites(ElectronDistributed const &block) { return block.nsites(); }
int64_t dim(ElectronDistributed const &block) { return block.dim(); }
int64_t size(ElectronDistributed const &block) { return block.size(); }
bool isreal(ElectronDistributed const &block) { return block.isreal(); }
std::ostream &operator<<(std::ostream &out, ElectronDistributed const &block) {
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  out << "ElectronDistributed:\n";
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
std::string to_string(ElectronDistributed const &block) {
  return to_string_generic(block);
}

ElectronDistributedIterator::ElectronDistributedIterator(
    ElectronDistributed const &block, bool begin)
    : nsites_(block.nsites()), pstate_(nsites_),
      it_(std::visit(
          [&](auto const &basis) {
            basis::BasisElectronDistributedIterator it =
                begin ? basis.begin() : basis.end();
            return it;
          },
          block.basis())) {}

ElectronDistributedIterator &ElectronDistributedIterator::operator++() {
  std::visit([](auto &&it) { ++it; }, it_);
  return *this;
}

ProductState const &ElectronDistributedIterator::operator*() const {
  std::visit(
      [&](auto &&it) {
        auto [ups, dns] = *it;
        to_product_state_electron(ups, dns, pstate_);
      },
      it_);
  return pstate_;
}

bool ElectronDistributedIterator::operator!=(
    ElectronDistributedIterator const &rhs) const {
  return it_ != rhs.it_;
}

} // namespace xdiag
