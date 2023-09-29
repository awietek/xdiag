#ifdef HYDRA_ENABLE_MPI
#include "basis_np.h"

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>

namespace hydra::basis::tj_distributed {

template <typename bit_t>
BasisNp::BasisNp<bit_t>(int64_t n_sites, int64_t n_up, int64_t n_dn)
    : n_sites_(n_sites), n_up_(n_up), n_dn_(n_dn) {
  using namespace combinatorics;

  if (n_sites < 0) {
    throw(std::invalid_argument("n_sites < 0"));
  } else if ((n_up < 0) || (n_dn < 0)) {
    throw(std::invalid_argument("nup < 0 or ndn < 0"));
  } else if ((n_up + n_dn) > n_sites) {
    throw(std::invalid_argument("nup + ndn > n_sites"));
  }

  size_ = binomial(n_sites, n_up) * binomial(n_sites - n_up, n_dn);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
  sitesmask_ = ((bit_t)1 << n_sites) - 1;

  // ////////////////////////////////////////////////////////////////
  // Ordering  ups / dns

  // Determine local ups
  for (auto ups : Combinations<bit_t>(n_sites, n_up)) {
    if (rank(ups) == mpi_rank_) {
      my_ups_.append(ups);
    }
  }

  // Determine corresponding dns
  my_dns_for_ups_.reserve(my_ups_.size());
  int64_t my_dns_size = my_ups_.size() * binomial(n_sites - n_up, n_dn);
  my_dns_for_ups_storage_.reserve(my_dns_size);
  for (auto ups : my_ups_) {
    bit_t not_ups = (~ups) & sitesmask_;

    // Create all dns for my ups configurations
    int64_t dns_start = my_dns_for_ups_storage_.size();
    for (auto dnsc : Combinations<bit_t>(n_sites - n_up, n_dn)) {
      bit_t dns = bits::deposit(dnsc, not_ups);
      my_dns_for_ups_storage_.push_back(dns);
    }
    int64_t dns_end = my_dns_for_ups_storage_.size();

    my_dns_for_ups_.push_back(
        span(my_dns_for_ups_storage_.data() + dns_start, dns_end - dns_start));
  }

  // compute forward transpose communicator going from ups / dns to dns / ups
  std::vector<int64_t> n_states_i_send(mpi_size_, 0);
  for (int64_t i = 0; i < my_ups_.size(); ++i) {
    bit_t ups = my_ups_[i];
    for (auto dns : my_dns_for_ups_[i]) {
      int target = rank(dns);
      ++n_states_i_send[target];
    }
  }
  transpose_communicator_ = mpi::Communicator(n_states_i_send);

  // ////////////////////////////////////////////////////////////////
  // Ordering  dns / ups

  // Determine local dns
  for (auto dns : Combinations<bit_t>(n_sites, n_dn)) {
    if (rank(dns) == mpi_rank_) {
      my_dns_.append(dns);
    }
  }

  // Determine corresponding ups
  my_ups_for_dns_.reserve(my_dns_.size());
  int64_t my_ups_size = my_dns_.size() * binomial(n_sites - n_dn, n_up);
  my_ups_for_dns_storage_.reserve(my_ups_size);
  for (auto dns : my_dns_) {
    bit_t not_dns = (~dns) & sitesmask_;

    // Create all ups for my dns configurations
    int64_t ups_start = my_ups_for_dns_storage_.size();
    for (auto upsc : Combinations<bit_t>(n_sites - n_dn, n_up)) {
      bit_t ups = bits::deposit(upsc, not_dns);
      my_ups_for_dns_storage_.push_back(ups);
    }
    int64_t ups_end = my_ups_for_dns_storage_.size();

    my_ups_for_dns_.push_back(
        span(my_ups_for_dns_storage_.data() + ups_start, ups_end - ups_start));
  }
}

#endif
