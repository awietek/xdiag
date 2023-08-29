#include "electron.h"

#include <hydra/utils/logger.h>

namespace hydra {

using namespace basis;

Electron::Electron(int64_t n_sites)
    : n_sites_(n_sites), charge_conserved_(false), charge_(undefined_qn),
      sz_conserved_(false), sz_(undefined_qn), n_up_(undefined_qn),
      n_dn_(undefined_qn), symmetric_(false), permutation_group_(), irrep_() {

  if (n_sites < 0) {
    throw(std::invalid_argument("n_sites < 0"));
  }
  try {
    if (n_sites < 16) {
      basis_ =
          std::make_shared<basis_t>(electron::BasisNoNp<uint16_t>(n_sites));
    } else if (n_sites < 32) {
      basis_ =
          std::make_shared<basis_t>(electron::BasisNoNp<uint32_t>(n_sites));
    } else if (n_sites < 64) {
      basis_ =
          std::make_shared<basis_t>(electron::BasisNoNp<uint64_t>(n_sites));
    } else {
      throw(std::runtime_error(
          "blocks with more than 64 sites currently not implemented"));
    }
    size_ = hydra::size(*basis_);
  } catch (...) {
    rethrow(__func__, "Cannot create Basis for Electron");
  }
}

Electron::Electron(int64_t n_sites, int64_t nup, int64_t ndn)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      symmetric_(false), permutation_group_(), irrep_() {

  if (n_sites < 0) {
    throw(std::invalid_argument("n_sites < 0"));
  } else if ((nup < 0) || (nup > n_sites)) {
    throw(std::invalid_argument("Invalid value of nup"));
  } else if ((ndn < 0) || (ndn > n_sites)) {
    throw(std::invalid_argument("Invalid value of ndn"));
  }

  try {
    if (n_sites < 16) {
      basis_ = std::make_shared<basis_t>(
          electron::BasisNp<uint16_t>(n_sites, nup, ndn));
    } else if (n_sites < 32) {
      basis_ = std::make_shared<basis_t>(
          electron::BasisNp<uint32_t>(n_sites, nup, ndn));
    } else if (n_sites < 64) {
      basis_ = std::make_shared<basis_t>(
          electron::BasisNp<uint64_t>(n_sites, nup, ndn));
    } else {
      throw(std::runtime_error(
          "blocks with more than 64 sites currently not implemented"));
    }
    size_ = hydra::size(*basis_);
  } catch (...) {
    rethrow(__func__, "Cannot create Basis for Electron");
  }
}

Electron::Electron(int64_t n_sites, PermutationGroup group,
                   Representation irrep)
    : n_sites_(n_sites), charge_conserved_(false), charge_(undefined_qn),
      sz_conserved_(false), sz_(undefined_qn), n_up_(undefined_qn),
      n_dn_(undefined_qn), symmetric_(true),
      permutation_group_(allowed_subgroup(group, irrep)), irrep_(irrep) {

  if (n_sites < 0) {
    throw(std::invalid_argument("n_sites < 0"));
  } else if (n_sites != group.n_sites()) {
    throw(std::logic_error(
        "n_sites does not match the n_sites in PermutationGroup"));
  } else if (permutation_group_.size() != irrep.size()) {
    throw(std::logic_error("PermutationGroup and Representation do not have "
                           "same number of elements"));
  }

  try {
    if (n_sites < 16) {
      basis_ = std::make_shared<basis_t>(
          electron::BasisSymmetricNoNp<uint16_t>(n_sites, group, irrep));
    } else if (n_sites < 32) {
      basis_ = std::make_shared<basis_t>(
          electron::BasisSymmetricNoNp<uint32_t>(n_sites, group, irrep));
    } else if (n_sites < 64) {
      basis_ = std::make_shared<basis_t>(
          electron::BasisSymmetricNoNp<uint64_t>(n_sites, group, irrep));
    } else {
      throw(std::runtime_error(
          "blocks with more than 64 sites currently not implemented"));
    }
    size_ = hydra::size(*basis_);
  } catch (...) {
    rethrow(__func__, "Cannot create Basis for Electron");
  }
}

Electron::Electron(int64_t n_sites, int64_t nup, int64_t ndn,
                   PermutationGroup group, Representation irrep)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      symmetric_(true), permutation_group_(allowed_subgroup(group, irrep)),
      irrep_(irrep) {

  if (n_sites < 0) {
    throw(std::invalid_argument("n_sites < 0"));
  } else if ((nup < 0) || (nup > n_sites)) {
    throw(std::invalid_argument("Invalid value of nup"));
  } else if ((ndn < 0) || (ndn > n_sites)) {
    throw(std::invalid_argument("Invalid value of ndn"));
  } else if (n_sites != group.n_sites()) {
    throw(std::logic_error(
        "n_sites does not match the n_sites in PermutationGroup"));
  } else if (permutation_group_.size() != irrep.size()) {
    throw(std::logic_error("PermutationGroup and Representation do not have "
                           "same number of elements"));
  }

  try {
    if (n_sites < 16) {
      basis_ = std::make_shared<basis_t>(electron::BasisSymmetricNp<uint16_t>(
          n_sites, nup, ndn, group, irrep));
    } else if (n_sites < 32) {
      basis_ = std::make_shared<basis_t>(electron::BasisSymmetricNp<uint32_t>(
          n_sites, nup, ndn, group, irrep));
    } else if (n_sites < 64) {
      basis_ = std::make_shared<basis_t>(electron::BasisSymmetricNp<uint64_t>(
          n_sites, nup, ndn, group, irrep));
    } else {
      throw(std::runtime_error(
          "blocks with more than 64 sites currently not implemented"));
    }
    size_ = hydra::size(*basis_);
  } catch (...) {
    rethrow(__func__, "Cannot create Basis for Electron");
  }
}

bool Electron::iscomplex(double precision) const {
  return symmetric_ ? irrep_.iscomplex(precision) : false;
}
bool Electron::isreal(double precision) const { return !iscomplex(precision); }

bool Electron::operator==(Electron const &rhs) const {
  return (n_sites_ == rhs.n_sites_) &&
         (charge_conserved_ == rhs.charge_conserved_) &&
         (charge_ == rhs.charge_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}
bool Electron::operator!=(Electron const &rhs) const {
  return !operator==(rhs);
}

basis_electron_variant_t const &Electron::basis() const { return *basis_; }

} // namespace hydra
