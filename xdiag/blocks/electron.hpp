#pragma once

#include <xdiag/basis/electron/basis_electron.hpp>
#include <xdiag/common.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

class Electron {
public:
  using basis_t = basis::BasisElectron;

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
  PermutationGroup const &permutation_group() const;
  Representation const &irrep() const;

  int64_t dim() const;
  int64_t size() const;
  bool isreal(double precision = 1e-12) const;

  bool operator==(Electron const &rhs) const;
  bool operator!=(Electron const &rhs) const;

  basis_t const &basis() const;

private:
  int64_t n_sites_;
  int64_t n_up_;
  int64_t n_dn_;
  PermutationGroup permutation_group_;
  Representation irrep_;
  std::shared_ptr<basis_t> basis_;
  int64_t size_;
};

} // namespace xdiag
