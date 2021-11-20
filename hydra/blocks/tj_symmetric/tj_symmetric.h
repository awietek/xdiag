#pragma once

#include <map>

#include <hydra/blocks/blocks.h>
#include <hydra/blocks/tj_symmetric/tj_symmetric_apply.h>
#include <hydra/blocks/tj_symmetric/tj_symmetric_indexing.h>
#include <hydra/blocks/tj_symmetric/tj_symmetric_matrix.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra {

template <typename bit_t> class tJSymmetric {
public:
  tJSymmetric() = default;
  tJSymmetric(int n_sites, int nup, int ndn, PermutationGroup permutation_group,
              Representation irrep);

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

  bool operator==(tJSymmetric<bit_t> const &rhs) const;
  bool operator!=(tJSymmetric<bit_t> const &rhs) const;

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
  int n_symmetries_;

  idx_t size_;

  indexing::tJSymmetricIndexing<bit_t> indexing_;
  indexing::tJSymmetricIndexing<bit_t> const &indexing() const {
    return indexing_;
  }
  
  friend lila::Matrix<double> MatrixReal<bit_t>(BondList const &,
                                                Couplings const &,
                                                tJSymmetric<bit_t> const &,
                                                tJSymmetric<bit_t> const &);
  friend lila::Matrix<complex> MatrixCplx<bit_t>(BondList const &,
                                                 Couplings const &,
                                                 tJSymmetric<bit_t> const &,
                                                 tJSymmetric<bit_t> const &);
  friend void Apply<bit_t>(BondList const &, Couplings const &,
                           tJSymmetric<bit_t> const &,
                           lila::Vector<double> const &,
                           tJSymmetric<bit_t> const &, lila::Vector<double> &);
  friend void Apply<bit_t>(BondList const &, Couplings const &,
                           tJSymmetric<bit_t> const &,
                           lila::Vector<complex> const &,
                           tJSymmetric<bit_t> const &, lila::Vector<complex> &);
};

} // namespace hydra
