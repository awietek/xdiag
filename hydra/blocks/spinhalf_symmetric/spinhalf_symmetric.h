#pragma once

#include <hydra/common.h>

#include <hydra/blocks/spinhalf_symmetric/spinhalf_symmetric_apply.h>
#include <hydra/blocks/spinhalf_symmetric/spinhalf_symmetric_matrix.h>

#include <hydra/indexing/spinhalf/spinhalf_symmetric_indexing.h>

#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <typename bit_t> class SpinhalfSymmetric {
public:
  SpinhalfSymmetric() = default;
  SpinhalfSymmetric(int n_sites, int n_up, PermutationGroup permutation_group,
                    Representation irrep);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_up_; }
  inline bool sz_conserved() const { return sz_conserved_; }

  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }
  inline Representation const &irrep() const { return irrep_; }
  inline idx_t size() const { return indexing_.size(); }

  bool operator==(SpinhalfSymmetric const &rhs) const;
  bool operator!=(SpinhalfSymmetric const &rhs) const;

  indexing::SpinhalfSymmetricIndexing<bit_t> const &indexing() const {
    return indexing_;
  }

private:
  int n_sites_;
  bool sz_conserved_;
  int n_up_;
  int n_dn_;
  int sz_;

  PermutationGroup permutation_group_;
  Representation irrep_;
  indexing::SpinhalfSymmetricIndexing<bit_t> indexing_;
};

} // namespace hydra
