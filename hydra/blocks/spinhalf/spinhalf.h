#pragma once

#include <memory>

#include <hydra/common.h>

#include <hydra/indexing/indexing_variants.h>

#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <typename bit_t> class Spinhalf {
public:
  Spinhalf() = default;
  Spinhalf(int n_sites);
  Spinhalf(int n_sites, int n_up);
  Spinhalf(int n_sites, PermutationGroup permutation_group,
           Representation irrep);
  Spinhalf(int n_sites, int n_up, PermutationGroup permutation_group,
           Representation irrep);

  inline int n_sites() const { return n_sites_; }
  inline bool sz_conserved() const { return sz_conserved_; }
  inline int sz() const { return sz_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }
  inline bool symmetric() const { return symmetric_; }
  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }
  inline Representation const &irrep() const { return irrep_; }
  inline idx_t size() const { return size_; }

  bool operator==(Spinhalf const &rhs) const;
  bool operator!=(Spinhalf const &rhs) const;

private:
  int n_sites_;
  bool sz_conserved_;
  int n_up_;
  int n_dn_;
  int sz_;
  bool symmetric_;
  PermutationGroup permutation_group_;
  Representation irrep_;

  using indexing_t = indexing::SpinhalfIndexing<bit_t>;
  std::shared_ptr<indexing_t> indexing_;
  indexing_t const &indexing() const;

  idx_t size_;

  template <typename bit_tt, typename coeff_tt>
  friend void
  Apply(BondList const &bonds, Couplings const &couplings,
        Spinhalf<bit_tt> const &block_in, lila::Vector<coeff_tt> const &vec_in,
        Spinhalf<bit_tt> const &block_out, lila::Vector<coeff_tt> &vec_out);

  template <typename bit_tt, typename coeff_tt>
  friend lila::Matrix<coeff_tt> MatrixGen(BondList const &bonds,
                                          Couplings const &couplings,
                                          Spinhalf<bit_tt> const &block_in,
                                          Spinhalf<bit_tt> const &block_out);
};

} // namespace hydra
