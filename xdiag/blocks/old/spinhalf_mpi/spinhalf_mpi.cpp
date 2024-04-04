#include "spinhalf_mpi.hpp"

#include <mpi.h>

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/hashes.hpp>
#include <xdiag/combinatorics/subsets.hpp>

#include <xdiag/mpi/allreduce.hpp>

#include <xdiag/bitops/bitops.hpp>
#include <xdiag/blocks/utils/block_utils.hpp>

namespace xdiag {

template <class bit_t>
SpinhalfMPI<bit_t>::SpinhalfMPI(int n_sites, int n_up)
    : n_sites_(n_sites), sz_conserved_(true), n_up_(n_up),
      n_dn_(n_sites - n_up), sz_(n_up_ - n_dn_),
      indexing_(std::make_shared<indexing_t>(indexing_t(n_sites, n_up))),
      size_(indexing_->size()), dim_(indexing_->dim()) {

  utils::check_nup_spinhalf(n_sites, n_up, "SpinhalfMPI");

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

  assert(dim_ == combinatorics::binomial(n_sites, n_up));
}

template <class bit_t>
bool SpinhalfMPI<bit_t>::operator==(SpinhalfMPI<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_);
}

template <class bit_t>
bool SpinhalfMPI<bit_t>::operator!=(SpinhalfMPI<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class SpinhalfMPI<uint16_t>;
template class SpinhalfMPI<uint32_t>;
template class SpinhalfMPI<uint64_t>;

} // namespace xdiag
