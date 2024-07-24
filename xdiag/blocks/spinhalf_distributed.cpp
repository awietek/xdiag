#include "spinhalf_distributed.hpp"

namespace xdiag {

SpinhalfDistributed::SpinhalfDistributed(int64_t n_sites, int64_t n_up) try
    : n_sites_(n_sites), n_up_(n_up) {
  using namespace basis::spinhalf_distributed;
  using combinatorics::binomial;

  if (n_sites < 0) {
    XDIAG_THROW("n_sites < 0");
  } else if (n_up < 0) {
    XDIAG_THROW("n_up < 0");
  } else if (n_up > n_sites) {
    XDIAG_THROW("n_up > n_sites");
  }

  if (n_sites < 32) {
    basis_ = std::make_shared<basis_t>(BasisSz<uint32_t>(n_sites, n_up));
  } else if (n_sites < 64) {
    basis_ = std::make_shared<basis_t>(BasisSz<uint64_t>(n_sites, n_up));
  } else {
    XDIAG_THROW("blocks with more than 64 sites currently not implemented");
  }
  dim_ = basis::dim(*basis_);
  assert(dim_ == binomial(n_sites, n_up));
  size_ = basis::size(*basis_);

  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t SpinhalfDistributed::n_sites() const { return n_sites_; }
int64_t SpinhalfDistributed::n_up() const { return n_up_; }

int64_t SpinhalfDistributed::dim() const { return dim_; }
int64_t SpinhalfDistributed::size() const { return size_; }
int64_t SpinhalfDistributed::size_max() const {
  return basis::size_max(*basis_);
}
int64_t SpinhalfDistributed::size_min() const {
  return basis::size_min(*basis_);
}

bool SpinhalfDistributed::isreal(double precision) const {
  return true; // would only be nontrivial with space group irreps
}

bool SpinhalfDistributed::operator==(SpinhalfDistributed const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_up_ == rhs.n_up_);
}
bool SpinhalfDistributed::operator!=(SpinhalfDistributed const &rhs) const {
  return !operator==(rhs);
}

SpinhalfDistributed::basis_t const &SpinhalfDistributed::basis() const {
  return *basis_;
}

} // namespace xdiag
