#include "electron_mpi.hpp"

#include <mpi.h>

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/mpi/allreduce.hpp>

namespace xdiag {

template <class bit_t>
ElectronMPI<bit_t>::ElectronMPI(int n_sites, int nup, int ndn)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      lintable_up_(n_sites, n_up_), lintable_dn_(n_sites, n_dn_),
      size_up_(combinatorics::binomial(n_sites, n_up_)),
      size_dn_(combinatorics::binomial(n_sites, n_dn_)),
      size_(size_up_ * size_dn_)

{
  using combinatorics::Combinations;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
  mpi::Allreduce(&size_, &dim_, 1, MPI_SUM, MPI_COMM_WORLD);

  // Determine the upspin configurations of this process
  int64_t offset = 0;
  for (auto ups : Combinations<bit_t>(n_sites, nup))
    if (process(ups) == mpi_rank_) {
      my_ups_.push_back(ups);
      my_ups_offset_[ups] = offset;
      offset += size_dn_;
    }

  // Determine the dnspin configurations of this process
  offset = 0;
  for (auto dns : Combinations<bit_t>(n_sites, ndn))
    if (process(dns) == mpi_rank_) {
      my_dns_.push_back(dns);
      my_dns_offset_[dns] = offset;
      offset += size_up_;
    }
}

template <class bit_t>
bool ElectronMPI<bit_t>::operator==(ElectronMPI<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_);
}

template <class bit_t>
bool ElectronMPI<bit_t>::operator!=(ElectronMPI<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class ElectronMPI<uint16_t>;
template class ElectronMPI<uint32_t>;
template class ElectronMPI<uint64_t>;

} // namespace xdiag
