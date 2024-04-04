#pragma once

#include <xdiag/basis/basis.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/bondlist.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

class Electron {
public:
  using basis_t = basis_electron_variant_t;

  Electron() = default;
  Electron(int64_t n_sites);
  Electron(int64_t n_sites, int64_t nup, int64_t ndn);
  Electron(int64_t n_sites, PermutationGroup permutation_group,
           Representation irrep);
  Electron(int64_t n_sites, int64_t nup, int64_t ndn,
           PermutationGroup permutation_group, Representation irrep);

  int64_t n_sites() const;
  int64_t n_up() const;
  int64_t n_dn() const;

  bool charge_conserved() const;
  bool sz_conserved() const;

  bool symmetric() const;
  PermutationGroup const &permutation_group() const;
  Representation const &irrep() const;

  int64_t dim() const;
  int64_t size() const;
  bool iscomplex(double precision = 1e-12) const;
  bool isreal(double precision = 1e-12) const;

  bool operator==(Electron const &rhs) const;
  bool operator!=(Electron const &rhs) const;

  basis_t const &basis() const;

private:
  int64_t n_sites_;
  bool charge_conserved_;
  int64_t charge_;
  bool sz_conserved_;
  int64_t sz_;
  int64_t n_up_;
  int64_t n_dn_;
  bool symmetric_;
  PermutationGroup permutation_group_;
  Representation irrep_;
  std::shared_ptr<basis_t> basis_;
  int64_t size_;
};

} // namespace xdiag
