#pragma once

#include <hydra/common.h>

#include <hydra/blocks/blocks.h>
#include <hydra/blocks/electron_symmetric/electron_symmetric_apply.h>
#include <hydra/blocks/electron_symmetric/electron_symmetric_matrix.h>

#include <hydra/indexing/electron/electron_symmetric_indexing.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra {

template <typename bit_t> class ElectronSymmetric {
public:
  ElectronSymmetric() = default;
  ElectronSymmetric(int n_sites, int nup, int ndn,
                    PermutationGroup permutation_group, Representation irrep);

  int n_sites() const { return n_sites_; }
  int n_up() const { return n_up_; }
  int n_dn() const { return n_up_; }
  bool charge_conserved() const { return charge_conserved_; }
  bool sz_conserved() const { return sz_conserved_; }

  PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }
  Representation const &irrep() const { return irrep_; }

  idx_t size() const { return size_; }

  bool operator==(ElectronSymmetric const &rhs) const;
  bool operator!=(ElectronSymmetric const &rhs) const;

  indexing::ElectronSymmetricIndexing<bit_t> const &indexing() const {
    return indexing_;
  }

private:
  int n_sites_;

  bool charge_conserved_;
  int charge_;
  bool sz_conserved_;
  int sz_;
  int n_up_;
  int n_dn_;

  PermutationGroup permutation_group_;
  Representation irrep_;

  indexing::ElectronSymmetricIndexing<bit_t> indexing_;
  idx_t size_;

};

} // namespace hydra
