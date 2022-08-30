#pragma once

#include <memory>

#include <lila/all.h>

#include <hydra/common.h>

#include <hydra/indexing/tj/tj_indexing.h>
#include <hydra/indexing/tj/tj_symmetric_indexing.h>

#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t = std_bit_t> class tJ {
public:
  tJ() = default;
  tJ(int n_sites, int nup, int ndn);
  tJ(int n_sites, int nup, int ndn, PermutationGroup permutation_group,
     Representation irrep);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }

  inline bool charge_conserved() const { return charge_conserved_; }
  inline bool sz_conserved() const { return sz_conserved_; }

  inline bool symmetric() const { return symmetric_; }
  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }
  inline Representation const &irrep() const { return irrep_; }
  inline idx_t size() const { return size_; }

  bool operator==(tJ const &rhs) const;
  bool operator!=(tJ const &rhs) const;

private:
  int n_sites_;

  bool charge_conserved_;
  int charge_;
  bool sz_conserved_;
  int sz_;
  int n_up_;
  int n_dn_;
  bool symmetric_;
  PermutationGroup permutation_group_;
  Representation irrep_;

  using indexing_np_t = indexing::tJIndexing<bit_t>;
  using indexing_sym_np_t = indexing::tJSymmetricIndexing<bit_t>;

  std::shared_ptr<indexing_np_t> indexing_np_;
  std::shared_ptr<indexing_sym_np_t> indexing_sym_np_;

  idx_t size_;

  indexing_np_t const &indexing_np() const;
  indexing_sym_np_t const &indexing_sym_np() const;

  template <typename bit_tt, typename coeff_tt>
  friend void
  Apply(BondList const &bonds, Couplings const &couplings,
        tJ<bit_tt> const &block_in, lila::Vector<coeff_tt> const &vec_in,
        tJ<bit_tt> const &block_out, lila::Vector<coeff_tt> &vec_out);

  template <typename bit_tt, typename coeff_tt>
  friend lila::Matrix<coeff_tt>
  MatrixGen(BondList const &bonds, Couplings const &couplings,
            tJ<bit_tt> const &block_in, tJ<bit_tt> const &block_out);
};

} // namespace hydra
