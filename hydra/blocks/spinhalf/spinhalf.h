#pragma once

#include <memory>

#include <hydra/common.h>

#include <hydra/blocks/spinhalf/spinhalf_apply.h>
#include <hydra/blocks/spinhalf/spinhalf_matrix.h>

#include <hydra/indexing/spinhalf/spinhalf_indexing.h>
#include <hydra/indexing/spinhalf/spinhalf_indexing_no_sz.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <typename bit_t> class Spinhalf {
public:
  Spinhalf() = default;
  Spinhalf(int n_sites);
  Spinhalf(int n_sites, int n_up);

  inline int n_sites() const { return n_sites_; }
  inline bool sz_conserved() const { return sz_conserved_; }
  inline int sz() const { return sz_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }
  inline idx_t size() const { return size_; }

  bool operator==(Spinhalf const &rhs) const;
  bool operator!=(Spinhalf const &rhs) const;

private:
  int n_sites_;
  bool sz_conserved_;
  int n_up_;
  int n_dn_;
  int sz_;
  idx_t size_;
  
  using idxng_sz_t = indexing::SpinhalfIndexing<bit_t>;
  using idxng_no_sz_t = indexing::SpinhalfIndexingNoSz<bit_t>;

  std::shared_ptr<idxng_sz_t> indexing_sz_conserved_;
  std::shared_ptr<idxng_no_sz_t> indexing_sz_not_conserved_;

  idxng_sz_t const &indexing_sz_conserved() const;
  idxng_no_sz_t const &indexing_sz_not_conserved() const;

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
