#pragma once
#include "extern/armadillo/armadillo"

#include <hydra/basis/basis.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra {

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

  inline int64_t n_sites() const { return n_sites_; }
  inline int64_t n_up() const { return n_up_; }
  inline int64_t n_dn() const { return n_dn_; }

  inline bool charge_conserved() const { return charge_conserved_; }
  inline bool sz_conserved() const { return sz_conserved_; }

  inline bool symmetric() const { return symmetric_; }
  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }
  inline Representation const &irrep() const { return irrep_; }

  inline int64_t size() const { return size_; }
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

} // namespace hydra
