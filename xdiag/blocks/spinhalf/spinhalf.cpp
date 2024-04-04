#include "spinhalf.hpp"

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag {

using namespace basis;

Spinhalf::Spinhalf(int64_t n_sites)
    : n_sites_(n_sites), sz_conserved_(false), n_up_(undefined_qn),
      n_dn_(undefined_qn), sz_(undefined_qn), symmetric_(false), n_sublat_(0),
      permutation_group_(), irrep_(), size_((int64_t)1 << n_sites) {
  try {

    if (n_sites < 0) {
      XDiagThrow(std::invalid_argument, "n_sites < 0");
    }

    if (n_sites < 32) {
      basis_ =
          std::make_shared<basis_t>(spinhalf::BasisNoSz<uint32_t>(n_sites));
    } else if (n_sites < 64) {
      basis_ =
          std::make_shared<basis_t>(spinhalf::BasisNoSz<uint64_t>(n_sites));
    } else {
      XDiagThrow(std::runtime_error,
                 "blocks with more than 64 sites currently not implemented");
    }

  } catch (...) {
    XDiagRethrow("Cannot create Basis for Spinhalf");
  }
}

Spinhalf::Spinhalf(int64_t n_sites, int64_t n_up)
    : n_sites_(n_sites), sz_conserved_(true), n_up_(n_up),
      n_dn_(n_sites - n_up), sz_(n_up_ - n_dn_), symmetric_(false),
      n_sublat_(0), permutation_group_(), irrep_(),
      size_(combinatorics::binomial(n_sites, n_up)) {
  try {

    if (n_sites < 0) {
      XDiagThrow(std::invalid_argument, "n_sites < 0");
    } else if ((n_up < 0) || (n_up > n_sites)) {
      XDiagThrow(std::invalid_argument, "Invalid value of nup");
    }

    if (n_sites < 32) {
      basis_ =
          std::make_shared<basis_t>(spinhalf::BasisSz<uint32_t>(n_sites, n_up));
    } else if (n_sites < 64) {
      basis_ =
          std::make_shared<basis_t>(spinhalf::BasisSz<uint64_t>(n_sites, n_up));
    } else {
      XDiagThrow(std::runtime_error,
                 "blocks with more than 64 sites currently not implemented");
    }

  } catch (...) {
    XDiagRethrow("Cannot create Basis for Spinhalf");
  }
}

template <typename bit_t>
std::shared_ptr<basis_spinhalf_variant_t>
make_spinhalf_basis_no_sz(int64_t n_sites, PermutationGroup const &group,
                          Representation const &irrep, int64_t n_sublat) try {

  std::shared_ptr<basis_spinhalf_variant_t> basis;
  if (n_sublat == 0) {
    basis = std::make_shared<basis_spinhalf_variant_t>(
        spinhalf::BasisSymmetricNoSz<bit_t>(n_sites, group, irrep));
  } else if (n_sublat == 1) {
    basis = std::make_shared<basis_spinhalf_variant_t>(
        spinhalf::BasisSublattice<bit_t, 1>(n_sites, group, irrep));
  } else if (n_sublat == 2) {
    basis = std::make_shared<basis_spinhalf_variant_t>(
        spinhalf::BasisSublattice<bit_t, 2>(n_sites, group, irrep));
  } else if (n_sublat == 3) {
    basis = std::make_shared<basis_spinhalf_variant_t>(
        spinhalf::BasisSublattice<bit_t, 3>(n_sites, group, irrep));
  } else if (n_sublat == 4) {
    basis = std::make_shared<basis_spinhalf_variant_t>(
        spinhalf::BasisSublattice<bit_t, 4>(n_sites, group, irrep));
  } else if (n_sublat == 5) {
    basis = std::make_shared<basis_spinhalf_variant_t>(
        spinhalf::BasisSublattice<bit_t, 5>(n_sites, group, irrep));
  } else {
    XDiagThrow(std::runtime_error,
               "Invalid n_sublat specified. Must be "
               "eiter 0 (so sublattice coding) or between 1 and 5. ");
  }
  return basis;
} catch (...) {
  XDiagRethrow("Cannot create Basis for Spinhalf without sz conservation but "
               "with permutation symmetries");
  return std::shared_ptr<basis_spinhalf_variant_t>();
}

Spinhalf::Spinhalf(int64_t n_sites, PermutationGroup group,
                   Representation irrep, int64_t n_sublat)
    : n_sites_(n_sites), sz_conserved_(false), n_up_(undefined_qn),
      n_dn_(undefined_qn), sz_(undefined_qn), symmetric_(true),
      n_sublat_(n_sublat), permutation_group_(allowed_subgroup(group, irrep)),
      irrep_(irrep) {
  try {

    if (n_sites < 0) {
      XDiagThrow(std::invalid_argument, "n_sites < 0");
    } else if (n_sites != group.n_sites()) {
      XDiagThrow(std::logic_error,
                 "n_sites does not match the n_sites in PermutationGroup");
    } else if (permutation_group_.size() != irrep.size()) {
      XDiagThrow(std::logic_error,
                 "PermutationGroup and Representation do not have "
                 "same number of elements");
    } else if ((n_sublat < 0) || (n_sublat > 5)) {
      XDiagThrow(std::invalid_argument,
                 "number of sublattices must either be 0 (no "
                 "sublattice) or between 1 and 5");
    }

    if (n_sites < 64) {
      basis_ =
          make_spinhalf_basis_no_sz<uint64_t>(n_sites, group, irrep, n_sublat);
    } else {
      throw(std::runtime_error(
          "blocks with more than 64 sites currently not implemented"));
    }
    size_ = xdiag::size(*basis_);

  } catch (...) {
    XDiagRethrow("Cannot create Basis for Spinhalf");
  }
}

template <typename bit_t>
std::shared_ptr<basis_spinhalf_variant_t>
make_spinhalf_basis_sz(int64_t n_sites, int64_t n_up,
                       PermutationGroup const &group,
                       Representation const &irrep, int64_t n_sublat) try {
  std::shared_ptr<basis_spinhalf_variant_t> basis;
  if (n_sublat == 0) {
    basis = std::make_shared<basis_spinhalf_variant_t>(
        spinhalf::BasisSymmetricSz<bit_t>(n_sites, n_up, group, irrep));
  } else if (n_sublat == 1) {
    basis = std::make_shared<basis_spinhalf_variant_t>(
        spinhalf::BasisSublattice<bit_t, 1>(n_sites, n_up, group, irrep));
  } else if (n_sublat == 2) {
    basis = std::make_shared<basis_spinhalf_variant_t>(
        spinhalf::BasisSublattice<bit_t, 2>(n_sites, n_up, group, irrep));
  } else if (n_sublat == 3) {
    basis = std::make_shared<basis_spinhalf_variant_t>(
        spinhalf::BasisSublattice<bit_t, 3>(n_sites, n_up, group, irrep));
  } else if (n_sublat == 4) {
    basis = std::make_shared<basis_spinhalf_variant_t>(
        spinhalf::BasisSublattice<bit_t, 4>(n_sites, n_up, group, irrep));
  } else if (n_sublat == 5) {
    basis = std::make_shared<basis_spinhalf_variant_t>(
        spinhalf::BasisSublattice<bit_t, 5>(n_sites, n_up, group, irrep));
  } else {
    XDiagThrow(std::invalid_argument,
               "Invalid n_sublat specified. Must be "
               "eiter 0 (so sublattice coding) or between 1 and 5.");
  }
  return basis;
} catch (...) {
  XDiagRethrow("Cannot create Basis for Spinhalf witt sz conservation and "
               "with permutation symmetries");
  return std::shared_ptr<basis_spinhalf_variant_t>();
}

Spinhalf::Spinhalf(int64_t n_sites, int64_t n_up, PermutationGroup group,
                   Representation irrep, int64_t n_sublat)
    : n_sites_(n_sites), sz_conserved_(true), n_up_(n_up),
      n_dn_(n_sites - n_up), sz_(n_up_ - n_dn_), symmetric_(true),
      n_sublat_(n_sublat), permutation_group_(allowed_subgroup(group, irrep)),
      irrep_(irrep) {
  try {
    if (n_sites < 0) {
      XDiagThrow(std::invalid_argument, "n_sites < 0");
    } else if ((n_up < 0) || (n_up > n_sites)) {
      XDiagThrow(std::invalid_argument, "Invalid value of nup");
    } else if (n_sites != group.n_sites()) {
      XDiagThrow(std::logic_error,
                 "n_sites does not match the n_sites in PermutationGroup");
    } else if (permutation_group_.size() != irrep.size()) {
      XDiagThrow(std::logic_error,
                 "PermutationGroup and Representation do not have "
                 "same number of elements");
    } else if ((n_sublat < 0) || (n_sublat > 5)) {
      XDiagThrow(std::invalid_argument,
                 "number of sublattices must either be 0 (no "
                 "sublattice) or between 1 and 5");
    }

    if (n_sites < 64) {
      basis_ = make_spinhalf_basis_sz<uint64_t>(n_sites, n_up, group, irrep,
                                                n_sublat);
    } else {
      XDiagThrow(std::runtime_error,
                 "blocks with more than 64 sites currently not implemented");
    }
    size_ = xdiag::size(*basis_);

  } catch (...) {
    XDiagRethrow("Cannot create Basis for Spinhalf");
  }
}

int64_t Spinhalf::n_sites() const { return n_sites_; }
bool Spinhalf::sz_conserved() const { return sz_conserved_; }
int64_t Spinhalf::sz() const { return sz_; }
int64_t Spinhalf::n_up() const { return n_up_; }
int64_t Spinhalf::n_dn() const { return n_dn_; }

bool Spinhalf::symmetric() const { return symmetric_; }
PermutationGroup Spinhalf::permutation_group() const {
  return permutation_group_;
}
Representation Spinhalf::irrep() const { return irrep_; }

int64_t Spinhalf::dim() const { return size_; }
int64_t Spinhalf::size() const { return size_; }

bool Spinhalf::iscomplex(double precision) const {
  return symmetric_ ? irrep_.iscomplex(precision) : false;
}
bool Spinhalf::isreal(double precision) const { return !iscomplex(precision); }

bool Spinhalf::operator==(Spinhalf const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (symmetric_ == rhs.symmetric_) && (n_sublat_ == rhs.n_sublat_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}

bool Spinhalf::operator!=(Spinhalf const &rhs) const {
  return !operator==(rhs);
}

basis_spinhalf_variant_t const &Spinhalf::basis() const { return *basis_; }

} // namespace xdiag
