#include "tj_distributed.hpp"

#include <xdiag/combinatorics/binomial.hpp>

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
} // namespace xdiag
