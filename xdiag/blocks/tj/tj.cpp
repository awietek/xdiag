#include "tj.h"

#include <xdiag/utils/logger.h>

namespace xdiag {

using namespace basis;

tJ::tJ(int64_t n_sites, int64_t nup, int64_t ndn)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      symmetric_(false), permutation_group_(), irrep_() {

  try {
    if (n_sites < 0) {
      XDiagThrow(std::invalid_argument, "n_sites < 0");
    } else if ((nup < 0) || (ndn < 0)) {
      XDiagThrow(std::invalid_argument, "nup < 0 or ndn < 0");
    } else if ((nup + ndn) > n_sites) {
      XDiagThrow(std::invalid_argument, "nup + ndn > n_sites");
    }

    if (n_sites < 16) {
      basis_ =
          std::make_shared<basis_t>(tj::BasisNp<uint16_t>(n_sites, nup, ndn));
    } else if (n_sites < 32) {
      basis_ =
          std::make_shared<basis_t>(tj::BasisNp<uint32_t>(n_sites, nup, ndn));
    } else if (n_sites < 64) {
      basis_ =
          std::make_shared<basis_t>(tj::BasisNp<uint64_t>(n_sites, nup, ndn));
    } else {
      XDiagThrow(std::runtime_error,
                 "blocks with more than 64 sites currently not implemented");
    }
    size_ = xdiag::size(*basis_);
  } catch (...) {
    XDiagRethrow("Cannot create Basis for tJ");
  }
}

tJ::tJ(int64_t n_sites, int64_t nup, int64_t ndn, PermutationGroup group,
       Representation irrep)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      symmetric_(true), permutation_group_(allowed_subgroup(group, irrep)),
      irrep_(irrep) {
  try {
    if (n_sites < 0) {
      XDiagThrow(std::invalid_argument, "n_sites < 0");
    } else if ((nup < 0) || (ndn < 0)) {
      XDiagThrow(std::invalid_argument, "nup < 0 or ndn < 0");
    } else if ((nup + ndn) > n_sites) {
      XDiagThrow(std::invalid_argument, "nup + ndn > n_sites");
    } else if (n_sites != group.n_sites()) {
      XDiagThrow(std::logic_error,
                 "n_sites does not match the n_sites in PermutationGroup");
    } else if (permutation_group_.size() != irrep.size()) {
      XDiagThrow(std::logic_error,
                 "PermutationGroup and Representation do not have "
                 "same number of elements");
    }

    if (n_sites < 16) {
      basis_ = std::make_shared<basis_t>(
          tj::BasisSymmetricNp<uint16_t>(n_sites, nup, ndn, group, irrep));
    } else if (n_sites < 32) {
      basis_ = std::make_shared<basis_t>(
          tj::BasisSymmetricNp<uint32_t>(n_sites, nup, ndn, group, irrep));
    } else if (n_sites < 64) {
      basis_ = std::make_shared<basis_t>(
          tj::BasisSymmetricNp<uint64_t>(n_sites, nup, ndn, group, irrep));
    } else {
      XDiagThrow(std::runtime_error,
                 "blocks with more than 64 sites currently not implemented");
    }
    size_ = xdiag::size(*basis_);
  } catch (...) {
    XDiagRethrow("Cannot create Basis for tJ");
  }
}

int64_t tJ::n_sites() const { return n_sites_; }
int64_t tJ::n_up() const { return n_up_; }
int64_t tJ::n_dn() const { return n_dn_; }

bool tJ::charge_conserved() const { return charge_conserved_; }
bool tJ::sz_conserved() const { return sz_conserved_; }

bool tJ::symmetric() const { return symmetric_; }
PermutationGroup const &tJ::permutation_group() const {
  return permutation_group_;
}
Representation const &tJ::irrep() const { return irrep_; }

int64_t tJ::dim() const { return size_; }
int64_t tJ::size() const { return size_; }

bool tJ::iscomplex(double precision) const {
  return symmetric_ ? irrep_.iscomplex(precision) : false;
}
bool tJ::isreal(double precision) const { return !iscomplex(precision); }

bool tJ::operator==(tJ const &rhs) const {
  return (n_sites_ == rhs.n_sites_) &&
         (charge_conserved_ == rhs.charge_conserved_) &&
         (charge_ == rhs.charge_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}
bool tJ::operator!=(tJ const &rhs) const { return !operator==(rhs); }

basis_tj_variant_t const &tJ::basis() const { return *basis_; }

} // namespace xdiag
