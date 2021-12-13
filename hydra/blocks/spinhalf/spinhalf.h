#pragma once

#include <memory>

#include <hydra/common.h>

#include <hydra/blocks/spinhalf/spinhalf_apply.h>
#include <hydra/blocks/spinhalf/spinhalf_matrix.h>

#include <hydra/indexing/spinhalf/spinhalf_indexing.h>
#include <hydra/indexing/spinhalf/spinhalf_indexing_no_sz.h>
#include <hydra/indexing/spinhalf/spinhalf_symmetric_indexing.h>
#include <hydra/indexing/spinhalf/spinhalf_symmetric_indexing_no_sz.h>

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

  using idxng_sz_t = indexing::SpinhalfIndexing<bit_t>;
  using idxng_no_sz_t = indexing::SpinhalfIndexingNoSz<bit_t>;
  using idxng_sym_sz_t = indexing::SpinhalfSymmetricIndexing<bit_t>;
  using idxng_sym_no_sz_t = indexing::SpinhalfSymmetricIndexingNoSz<bit_t>;

  std::shared_ptr<idxng_sz_t> indexing_sz_conserved_;
  std::shared_ptr<idxng_no_sz_t> indexing_sz_not_conserved_;
  std::shared_ptr<idxng_sym_sz_t> indexing_sym_sz_conserved_;
  std::shared_ptr<idxng_sym_no_sz_t> indexing_sym_sz_not_conserved_;

  idx_t size_;

  idxng_sz_t const &indexing_sz_conserved() const;
  idxng_no_sz_t const &indexing_sz_not_conserved() const;
  idxng_sym_sz_t const &indexing_sym_sz_conserved() const;
  idxng_sym_no_sz_t const &indexing_sym_sz_not_conserved() const;

  friend void Apply<bit_t, double>(BondList const &bonds,
                                   Couplings const &couplings,
                                   Spinhalf<bit_t> const &block_in,
                                   lila::Vector<double> const &vec_in,
                                   Spinhalf<bit_t> const &block_out,
                                   lila::Vector<double> &vec_out);
  friend void Apply<bit_t, complex>(BondList const &bonds,
                                    Couplings const &couplings,
                                    Spinhalf<bit_t> const &block_in,
                                    lila::Vector<complex> const &vec_in,
                                    Spinhalf<bit_t> const &block_out,
                                    lila::Vector<complex> &vec_out);
  friend lila::Matrix<double>
  MatrixGen<bit_t, double>(BondList const &bonds, Couplings const &couplings,
                           Spinhalf<bit_t> const &block_in,
                           Spinhalf<bit_t> const &block_out);
  friend lila::Matrix<complex>
  MatrixGen<bit_t, complex>(BondList const &bonds, Couplings const &couplings,
                            Spinhalf<bit_t> const &block_in,
                            Spinhalf<bit_t> const &block_out);
};

} // namespace hydra
