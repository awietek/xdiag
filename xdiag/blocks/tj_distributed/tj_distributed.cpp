#include "tj_distributed.hpp"

#include <xdiag/combinatorics/binomial.hpp>

namespace xdiag {

tJDistributed::tJDistributed(int64_t n_sites, int64_t n_up, int64_t n_dn)
    : n_sites_(n_sites), n_up_(n_up), n_dn_(n_dn) {
  using namespace basis::tj_distributed;
  using combinatorics::binomial;

  try {
    if (n_sites < 0) {
      XDIAG_THROW("n_sites < 0");
    } else if ((n_up < 0) || (n_dn < 0)) {
      XDIAG_THROW("n_up < 0 or n_dn < 0");
    } else if ((n_up + n_dn) > n_sites) {
      XDIAG_THROW("n_up + n_dn > n_sites");
    }

    if (n_sites < 16) {
      basis_ =
          std::make_shared<basis_t>(BasisNp<uint16_t>(n_sites, n_up, n_dn));
    } else if (n_sites < 32) {
      basis_ =
          std::make_shared<basis_t>(BasisNp<uint32_t>(n_sites, n_up, n_dn));
    } else if (n_sites < 64) {
      basis_ =
          std::make_shared<basis_t>(BasisNp<uint64_t>(n_sites, n_up, n_dn));
    } else {
      XDIAG_THROW("blocks with more than 64 sites currently not implemented");
    }
    dim_ = xdiag::dim(*basis_);
    assert(dim_ == binomial(n_sites, n_up) * binomial(n_sites - n_up, n_dn));
    size_ = xdiag::size(*basis_);
  } catch (Error const &e) {
    XDIAG_RETHROW(e);
  }
}

int64_t tJDistributed::n_sites() const { return n_sites_; }
int64_t tJDistributed::n_up() const { return n_up_; }
int64_t tJDistributed::n_dn() const { return n_dn_; }

int64_t tJDistributed::dim() const { return dim_; }
int64_t tJDistributed::size() const { return size_; }
int64_t tJDistributed::size_max() const { return xdiag::size_max(*basis_); }
int64_t tJDistributed::size_min() const { return xdiag::size_min(*basis_); }

bool tJDistributed::charge_conserved() const { return true; }
bool tJDistributed::sz_conserved() const { return true; }

bool tJDistributed::symmetric() const { return false; }
PermutationGroup tJDistributed::permutation_group() const {
  return PermutationGroup();
}
Representation tJDistributed::irrep() const { return Representation(); }

bool tJDistributed::iscomplex(double precision) const {
  (void)precision;
  return false;
}
bool tJDistributed::isreal(double precision) const {
  return !iscomplex(precision);
}

bool tJDistributed::operator==(tJDistributed const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_up_ == rhs.n_up_) &&
         (n_dn_ == rhs.n_dn_);
}
bool tJDistributed::operator!=(tJDistributed const &rhs) const {
  return !operator==(rhs);
}

basis_tj_distributed_variant_t const &tJDistributed::basis() const {
  return *basis_;
}
} // namespace xdiag
