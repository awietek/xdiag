#include "tj.hpp"

#include <xdiag/utils/logger.hpp>

namespace xdiag {

using namespace basis;

tJ::tJ(int64_t n_sites, int64_t nup, int64_t ndn) try
    : n_sites_(n_sites), n_up_(nup), n_dn_(ndn), permutation_group_(),
      irrep_() {

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
  size_ = xdiag::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

tJ::tJ(int64_t n_sites, int64_t nup, int64_t ndn, PermutationGroup group,
       Representation irrep) try
    : n_sites_(n_sites), n_up_(nup), n_dn_(ndn),
      permutation_group_(allowed_subgroup(group, irrep)), irrep_(irrep) {

  if (n_sites < 0) {
    XDIAG_THROW("Invalid argument: n_sites < 0");
  } else if ((nup < 0) || (ndn < 0)) {
    XDIAG_THROW("Invalid argument: (nup < 0) or (ndn < 0)");
  } else if ((nup + ndn) > n_sites) {
    XDIAG_THROW("Invalid argument: nup + ndn > n_sites");
  } else if (n_sites != group.n_sites()) {
    XDIAG_THROW("n_sites does not match the n_sites in PermutationGroup");
  } else if (permutation_group_.size() != irrep.size()) {
    XDIAG_THROW("PermutationGroup and Representation do not have "
                "same number of elements");
  } else if (n_sites < 32) {
    basis_ = std::make_shared<basis_t>(
        tj::BasisSymmetricNp<uint32_t>(n_sites, nup, ndn, group, irrep));
  } else if (n_sites < 64) {
    basis_ = std::make_shared<basis_t>(
        tj::BasisSymmetricNp<uint64_t>(n_sites, nup, ndn, group, irrep));
  } else {
    XDIAG_THROW("blocks with more than 64 sites currently not implemented");
  }
  size_ = xdiag::size(*basis_);
  check_dimension_works_with_blas_int_size(size_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t tJ::n_sites() const { return n_sites_; }
int64_t tJ::n_up() const { return n_up_; }
int64_t tJ::n_dn() const { return n_dn_; }

PermutationGroup const &tJ::permutation_group() const {
  return permutation_group_;
}
Representation const &tJ::irrep() const { return irrep_; }

int64_t tJ::dim() const { return size_; }
int64_t tJ::size() const { return size_; }

bool tJ::isreal(double precision) const { return irrep_.isreal(precision); }

bool tJ::operator==(tJ const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_up_ == rhs.n_up_) &&
         (n_dn_ == rhs.n_dn_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}
bool tJ::operator!=(tJ const &rhs) const { return !operator==(rhs); }

tJ::basis_t const &tJ::basis() const { return *basis_; }

} // namespace xdiag
