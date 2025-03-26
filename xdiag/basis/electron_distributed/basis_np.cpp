#include "basis_np.hpp"

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/parallel/mpi/allreduce.hpp>

namespace xdiag::basis::electron_distributed {

template <typename bit_t>
BasisNp<bit_t>::BasisNp(int64_t nsites, int64_t nup, int64_t ndn) try
    : nsites_(nsites), nup_(nup), ndn_(ndn), lintable_dns_(nsites, ndn),
      lintable_ups_(nsites, nup) {
  check_nsites_work_with_bits<bit_t>(nsites_);

  using namespace combinatorics;

  // Safety checks
  if (nsites < 0) {
    XDIAG_THROW("nsites < 0");
  } else if ((nup < 0) || (ndn < 0)) {
    XDIAG_THROW("nup < 0 or ndn < 0");
  } else if (nup > nsites) {
    XDIAG_THROW("nup > nsites");
  } else if (ndn > nsites) {
    XDIAG_THROW("ndn > nsites");
  }

  dim_ = binomial(nsites, nup) * binomial(nsites, ndn);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
  sitesmask_ = ((bit_t)1 << nsites) - 1;

  // ////////////////////////////////////////////////////////////////
  // Ordering  ups / dns

  // Determine local ups
  for (auto ups : Combinations<bit_t>(nsites, nup)) {
    if (rank(ups) == mpi_rank_) {
      my_ups_offset_[ups] = my_ups_.size() * binomial(nsites, ndn);
      my_ups_.push_back(ups);
    }
  }

  all_dns_.reserve(binomial(nsites, ndn));
  for (auto dns : Combinations<bit_t>(nsites, ndn)) {
    all_dns_.push_back(dns);
  }

  size_ = my_ups_.size() * all_dns_.size();

  // Determine local dns
  for (auto dns : Combinations<bit_t>(nsites, ndn)) {
    if (rank(dns) == mpi_rank_) {
      my_dns_offset_[dns] = my_dns_.size() * binomial(nsites, nup);
      my_dns_.push_back(dns);
    }
  }

  all_ups_.reserve(binomial(nsites, nup));
  for (auto ups : Combinations<bit_t>(nsites, nup)) {
    all_ups_.push_back(ups);
  }

  size_transpose_ = my_dns_.size() * all_ups_.size();

  ////
  // Determine transpose communicator patterns (forward)
  ///

  // compute forward transpose communicator going from ups / dns to dns / ups
  std::vector<int64_t> n_states_i_send(mpi_size_, 0);
  for (bit_t up : my_ups_) {
    for (bit_t dn : all_dns_) {
      int target = rank(dn);
      ++n_states_i_send[target];
    }
  }
  transpose_communicator_ = mpi::Communicator(n_states_i_send);
  int64_t send_size = transpose_communicator_.send_buffer_size();
  int64_t recv_size = transpose_communicator_.recv_buffer_size();
  mpi::buffer.reserve<bit_t>(send_size, recv_size);

  bit_t *send_buffer = mpi::buffer.send<bit_t>();

  // Send the up configurations
  transpose_communicator_.flush();
  std::vector<bit_t> ups_i_recv_by_transpose(recv_size, 0);
  for (bit_t up : my_ups_) {
    for (bit_t dn : all_dns_) {
      int target = rank(dn);
      transpose_communicator_.add_to_send_buffer(target, up, send_buffer);
    }
  }
  transpose_communicator_.all_to_all(send_buffer,
                                     ups_i_recv_by_transpose.data());
  // Send the dn configurations
  transpose_communicator_.flush();
  std::vector<bit_t> dns_i_recv_by_transpose(recv_size, 0);
  for (bit_t up : my_ups_) {
    for (bit_t dn : all_dns_) {
      int target = rank(dn);
      transpose_communicator_.add_to_send_buffer(target, dn, send_buffer);
    }
  }
  transpose_communicator_.all_to_all(send_buffer,
                                     dns_i_recv_by_transpose.data());

  // Compute the transpose permutation
  transpose_permutation_.resize(recv_size);
  for (int64_t i = 0; i < recv_size; ++i) {
    bit_t up = ups_i_recv_by_transpose[i];
    bit_t dn = dns_i_recv_by_transpose[i];
    int64_t idx = my_dns_offset_[dn];
    idx += index_ups(up);
    transpose_permutation_[i] = idx;
  }

  ////
  // Determine transpose communicator patterns (back)
  ////

  // compute forward transpose communicator going from ups / dns to dns / ups
  std::fill(n_states_i_send.begin(), n_states_i_send.end(), 0);
  for (bit_t dn : my_dns_) {
    for (bit_t up : all_ups_) {
      int target = rank(up);
      ++n_states_i_send[target];
    }
  }
  transpose_communicator_r_ = mpi::Communicator(n_states_i_send);
  send_size = transpose_communicator_r_.send_buffer_size();
  recv_size = transpose_communicator_r_.recv_buffer_size();
  mpi::buffer.reserve<bit_t>(send_size, recv_size);

  send_buffer = mpi::buffer.send<bit_t>();

  // Send the up configurations
  transpose_communicator_r_.flush();
  std::vector<bit_t> dns_i_recv_by_transpose_r(recv_size, 0);
  for (bit_t dn : my_dns_) {
    for (bit_t up : all_ups_) {
      int target = rank(up);
      transpose_communicator_r_.add_to_send_buffer(target, dn, send_buffer);
    }
  }
  transpose_communicator_r_.all_to_all(send_buffer,
                                       dns_i_recv_by_transpose_r.data());
  // Send the dn configurations
  transpose_communicator_r_.flush();
  std::vector<bit_t> ups_i_recv_by_transpose_r(recv_size, 0);
  for (bit_t dn : my_dns_) {
    for (bit_t up : all_ups_) {
      int target = rank(up);
      transpose_communicator_r_.add_to_send_buffer(target, up, send_buffer);
    }
  }
  transpose_communicator_r_.all_to_all(send_buffer,
                                       ups_i_recv_by_transpose_r.data());

  // Compute the transpose permutation
  transpose_permutation_r_.resize(recv_size);
  for (int64_t i = 0; i < recv_size; ++i) {
    bit_t up = ups_i_recv_by_transpose_r[i];
    bit_t dn = dns_i_recv_by_transpose_r[i];

    int64_t idx = my_ups_offset_[up];
    idx += index_dns(dn);
    transpose_permutation_r_[i] = idx;
  }

  transpose_communicator_.flush();
  transpose_communicator_r_.flush();

  // compute maximal size and size_transpose between processes
  int64_t size_max_f = 0;
  int64_t size_max_r = 0;
  mpi::Allreduce(&size_, &size_max_f, 1, MPI_MAX, MPI_COMM_WORLD);
  mpi::Allreduce(&size_transpose_, &size_max_r, 1, MPI_MAX, MPI_COMM_WORLD);
  size_max_ = std::max(size_max_f, size_max_r);

  int64_t size_min_f = 0;
  int64_t size_min_r = 0;
  mpi::Allreduce(&size_, &size_min_f, 1, MPI_MIN, MPI_COMM_WORLD);
  mpi::Allreduce(&size_transpose_, &size_min_r, 1, MPI_MIN, MPI_COMM_WORLD);
  size_min_ = std::min(size_min_f, size_min_r);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename bit_t> int64_t BasisNp<bit_t>::nsites() const {
  return nsites_;
}
template <typename bit_t> int64_t BasisNp<bit_t>::nup() const { return nup_; }

template <typename bit_t> int64_t BasisNp<bit_t>::ndn() const { return ndn_; }

template <typename bit_t>
int64_t BasisNp<bit_t>::index(bit_t up, bit_t dn) const {
  if (rank(up) != mpi_rank_) {
    return invalid_index;
  } else {
    return my_ups_offset(up) + index_dns(dn);
  }
}

template <typename bit_t> int64_t BasisNp<bit_t>::dim() const { return dim_; }
template <typename bit_t> int64_t BasisNp<bit_t>::size() const { return size_; }
template <typename bit_t> int64_t BasisNp<bit_t>::size_transpose() const {
  return size_transpose_;
}
template <typename bit_t> int64_t BasisNp<bit_t>::size_max() const {
  return size_max_;
}
template <typename bit_t> int64_t BasisNp<bit_t>::size_min() const {
  return size_min_;
}
template <typename bit_t>
typename BasisNp<bit_t>::iterator_t BasisNp<bit_t>::begin() const {
  return iterator_t(*this, true);
}
template <typename bit_t>
typename BasisNp<bit_t>::iterator_t BasisNp<bit_t>::end() const {
  return iterator_t(*this, false);
}

template <typename bit_t>
std::vector<bit_t> const &BasisNp<bit_t>::my_ups() const {
  return my_ups_;
}

template <typename bit_t>
int64_t BasisNp<bit_t>::my_ups_offset(bit_t ups) const {
  return my_ups_offset_.at(ups);
}

template <typename bit_t>
std::vector<bit_t> const &BasisNp<bit_t>::all_dns() const {
  return all_dns_;
}

template <typename bit_t>
std::vector<bit_t> const &BasisNp<bit_t>::my_dns() const {
  return my_dns_;
}
template <typename bit_t>
int64_t BasisNp<bit_t>::my_dns_offset(bit_t dns) const {
  return my_dns_offset_.at(dns);
}
template <typename bit_t>
std::vector<bit_t> const &BasisNp<bit_t>::all_ups() const {
  return all_ups_;
}

template <typename bit_t>
template <typename coeff_t>
void BasisNp<bit_t>::transpose(const coeff_t *in_vec, coeff_t *out_vec) const
    try {
  // transforms a vector in up/dn order to dn/up order
  // result of transpose is stored in send_buffer

  auto comm = transpose_communicator_;
  coeff_t *send_buffer = mpi::buffer.send<coeff_t>();
  coeff_t *recv_buffer = mpi::buffer.recv<coeff_t>();

  // Send the up dn configurations around
  int64_t idx = 0;
  for (bit_t up : my_ups_) {
    for (bit_t dn : all_dns_) {
      int target = rank(dn);
      comm.add_to_send_buffer(target, in_vec[idx], send_buffer);
      ++idx;
    }
  }

  assert(idx == size_);
  comm.all_to_all(send_buffer, recv_buffer);

  // Sort to proper order
  if (out_vec) {
    for (idx = 0; idx < size_transpose_; ++idx) {
      int64_t sorted_idx = transpose_permutation_[idx];
      out_vec[sorted_idx] = recv_buffer[idx];
    }
  } else {
    for (idx = 0; idx < size_transpose_; ++idx) {
      int64_t sorted_idx = transpose_permutation_[idx];
      send_buffer[sorted_idx] = recv_buffer[idx];
    }
  }
  mpi::buffer.clean_recv();
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void BasisNp<uint16_t>::transpose(const double *, double *) const;
template void BasisNp<uint16_t>::transpose(const complex *, complex *) const;
template void BasisNp<uint32_t>::transpose(const double *, double *) const;
template void BasisNp<uint32_t>::transpose(const complex *, complex *) const;
template void BasisNp<uint64_t>::transpose(const double *, double *) const;
template void BasisNp<uint64_t>::transpose(const complex *, complex *) const;

template <typename bit_t>
template <typename coeff_t>
void BasisNp<bit_t>::transpose_r(coeff_t const *in_vec, coeff_t *out_vec) const
    try {
  // transforms a vector in up/dn order to dn/up order
  // result of transpose is stored in send_buffer
  auto comm = transpose_communicator_r_;
  coeff_t *send_buffer = mpi::buffer.send<coeff_t>();
  coeff_t *recv_buffer = mpi::buffer.recv<coeff_t>();

  // Send the up dn configurations around
  int64_t idx = 0;
  for (bit_t dns : my_dns_) {
    for (bit_t up : all_ups_) {
      int target = rank(up);
      comm.add_to_send_buffer(target, in_vec[idx], send_buffer);
      ++idx;
    }
  }
  assert(idx == size_transpose_);

  comm.all_to_all(send_buffer, recv_buffer);

  // Sort to proper order
  if (out_vec) {
    for (idx = 0; idx < size_; ++idx) {
      int64_t sorted_idx = transpose_permutation_r_[idx];
      out_vec[sorted_idx] = recv_buffer[idx];
    }
  } else {
    for (idx = 0; idx < size_; ++idx) {
      int64_t sorted_idx = transpose_permutation_r_[idx];
      send_buffer[sorted_idx] = recv_buffer[idx];
    }
  }

  mpi::buffer.clean_recv();
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void BasisNp<uint16_t>::transpose_r(const double *, double *) const;
template void BasisNp<uint16_t>::transpose_r(const complex *, complex *) const;
template void BasisNp<uint32_t>::transpose_r(const double *, double *) const;
template void BasisNp<uint32_t>::transpose_r(const complex *, complex *) const;
template void BasisNp<uint64_t>::transpose_r(const double *, double *) const;
template void BasisNp<uint64_t>::transpose_r(const complex *, complex *) const;

template class BasisNp<uint32_t>;
template class BasisNp<uint64_t>;

template <typename bit_t>
BasisNpIterator<bit_t>::BasisNpIterator(BasisNp<bit_t> const &basis, bool begin)
    : basis_(basis), up_idx_(begin ? 0 : basis_.my_ups().size()), dn_idx_(0) {}

template <typename bit_t>
BasisNpIterator<bit_t> &BasisNpIterator<bit_t>::operator++() {
  ++dn_idx_;
  if (dn_idx_ == basis_.all_dns().size()) {
    dn_idx_ = 0;
    ++up_idx_;
    if (up_idx_ == basis_.my_ups().size()) {
      return *this;
    }
  }
  return *this;
}

template <typename bit_t>
std::pair<bit_t, bit_t> BasisNpIterator<bit_t>::operator*() const {
  bit_t ups = basis_.my_ups()[up_idx_];
  bit_t dns = basis_.all_dns()[dn_idx_];
  return {ups, dns};
}

template <typename bit_t>
bool BasisNpIterator<bit_t>::operator!=(
    BasisNpIterator<bit_t> const &rhs) const {
  return (up_idx_ != rhs.up_idx_) || (dn_idx_ != rhs.dn_idx_);
}

template class BasisNpIterator<uint32_t>;
template class BasisNpIterator<uint64_t>;

} // namespace xdiag::basis::electron_distributed
