#pragma once

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/fermi_sign.h>
#include <hydra/symmetries/representation.h>
#include <hydra/symmetries/spacegroup.h>

namespace hydra {

template <class bit_t, class SymmetryGroup = SpaceGroup<bit_t>>
class ElectronSymmetric {
public:
  ElectronSymmetric() = default;
  ElectronSymmetric(int n_sites, int charge, int sz,
                    SymmetryGroup symmetry_group, Representation irrep);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_up_; }
  inline bool charge_conserved() const { return charge_conserved_; }
  inline bool sz_conserved() const { return sz_conserved_; }

  inline SymmetryGroup symmetry_group() const { return symmetry_group_; }
  inline Representation irrep() const { return irrep_; }

  inline idx_t size() const { return size_; }
  inline double ups(idx_t idx) const { return ups_[idx]; }
  inline double dns(idx_t idx) const { return dns_[idx]; }
  inline double norm(idx_t idx) const { return norms_[idx]; }

  std::tuple<bit_t, bit_t> representative(bit_t ups, bit_t dns);

  bool operator==(ElectronSymmetric const &rhs);
  bool operator!=(ElectronSymmetric const &rhs);

private:
  int n_sites_;

  bool charge_conserved_;
  int charge_;
  bool sz_conserved_;
  int sz_;
  int nup_;
  int ndn_;

  SymmetryGroup symmetry_group_;
  Representation irrep_;
  std::vector<bit_t> ups_;
  std::vector<bit_t> dns_;
  std::vector<double> norms_;

  idx_t size_;
};

} // namespace hydra
