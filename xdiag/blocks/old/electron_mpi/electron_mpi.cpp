#include "electron_mpi.hpp"

#include <mpi.h>

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/mpi/allreduce.hpp>

namespace xdiag {

template <class bit_t>
ElectronMPI<bit_t>::ElectronMPI(int nsites, int nup, int ndn)
    : nsites_(nsites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), nup_(nup), ndn_(ndn),
      lintable_up_(nsites, nup_), lintable_dn_(nsites, ndn_),
      size_up_(combinatorics::binomial(nsites, nup_)),
      size_dn_(combinatorics::binomial(nsites, ndn_)),
      size_(size_up_ * size_dn_)

{
  using combinatorics::Combinations;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
  mpi::Allreduce(&size_, &dim_, 1, MPI_SUM, MPI_COMM_WORLD);

  // Determine the upspin configurations of this process
  int64_t offset = 0;
  for (auto ups : Combinations<bit_t>(nsites, nup))
    if (process(ups) == mpi_rank_) {
      my_ups_.push_back(ups);
      my_ups_offset_[ups] = offset;
      offset += size_dn_;
    }

  // Determine the dnspin configurations of this process
  offset = 0;
  for (auto dns : Combinations<bit_t>(nsites, ndn))
    if (process(dns) == mpi_rank_) {
      my_dns_.push_back(dns);
      my_dns_offset_[dns] = offset;
      offset += size_up_;
    }
}

template <class bit_t>
bool ElectronMPI<bit_t>::operator==(ElectronMPI<bit_t> const &rhs) const {
  return (nsites_ == rhs.nsites_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (nup_ == rhs.nup_) && (ndn_ == rhs.ndn_);
}

template <class bit_t>
bool ElectronMPI<bit_t>::operator!=(ElectronMPI<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class ElectronMPI<uint16_t>;
template class ElectronMPI<uint32_t>;
template class ElectronMPI<uint64_t>;

} // namespace xdiag
