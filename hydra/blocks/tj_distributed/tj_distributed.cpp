#include "tj_distributed.h"

namespace hydra {

tJDistributed::tJDistributed(int64_t n_sites, int64_t n_up, int64_t n_dn)
    : n_sites_(n_sites), n_up_(n_up), n_dn_(n_dn) {}

int64_t tJDistributed::n_sites() const { return n_sites_; }
int64_t tJDistributed::n_up() const { return n_up_; }
int64_t tJDistributed::n_dn() const { return n_dn_; }

int64_t tJDistributed::size() const { return size_; }
int64_t tJDistributed::local_size() const { return local_size_; }

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
} // namespace hydra
