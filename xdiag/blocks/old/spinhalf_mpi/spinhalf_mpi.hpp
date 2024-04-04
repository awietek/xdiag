#pragma once
#ifdef XDIAG_ENABLE_MPI

#include <mpi.h>
#include <xdiag/common.hpp>

#include <xdiag/indexing/spinhalf_mpi/spinhalf_mpi_indexing_sz.hpp>

#include <xdiag/operators/bondlist.hpp>
#include <xdiag/operators/couplings.hpp>

namespace xdiag {

template <class bit_t> class SpinhalfMPI {
public:
  SpinhalfMPI() = default;
  SpinhalfMPI(int n_sites, int n_up);

  inline int n_sites() const { return n_sites_; }

  inline bool sz_conserved() const { return sz_conserved_; }
  inline int sz() const { return sz_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }

  inline int64_t size() const { return size_; }
  inline int64_t dim() const { return dim_; }

  inline int mpi_rank() const { return mpi_rank_; }
  inline int mpi_size() const { return mpi_size_; }

  bool operator==(SpinhalfMPI const &rhs) const;
  bool operator!=(SpinhalfMPI const &rhs) const;

private:
  int n_sites_;
  bool sz_conserved_;
  int n_up_;
  int n_dn_;
  int sz_;

  using indexing_t = indexing::SpinhalfMPIIndexingSz<bit_t>;
  std::shared_ptr<indexing_t> indexing_;
  inline indexing_t const &indexing() const {return *indexing_; }

  int64_t size_;
  int64_t dim_;

  int mpi_rank_;
  int mpi_size_;

  template <typename bit_tt, typename coeff_tt>
  friend void Apply(BondList const &bonds, Couplings const &couplings,
                    SpinhalfMPI<bit_tt> const &block_in,
                    lila::Vector<coeff_tt> const &vec_in,
                    SpinhalfMPI<bit_tt> const &block_out,
                    lila::Vector<coeff_tt> &vec_out);
};

} // namespace xdiag
#endif
