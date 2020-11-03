#include "communicator.h"

#include <numeric>

namespace hydra {
namespace mpi {

template <class idx_t>
Communicator<idx_t>::Communicator(std::vector<idx_t> const &n_values_i_send)
    : n_values_prepared_(n_values_i_send.size(), 0),
      n_values_i_recv_(n_values_i_send.size(), 0),
      n_values_i_send_offsets_(n_values_i_send.size(), 0),
      n_values_i_recv_offsets_(n_values_i_send.size(), 0) {
  for (auto num : n_values_i_send)
    n_values_i_send_.push_back((int)num);

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
  assert(n_values_i_send.size() == mpi_size_);
  MPI_Alltoall(n_values_i_send_.data(), 1, MPI_INT,
               n_values_i_recv_.data(), 1, MPI_INT, MPI_COMM_WORLD);

  for (int i = 0; i < mpi_size_; ++i) {
    n_values_i_send_offsets_[i] = std::accumulate(
        n_values_i_send_.begin(), n_values_i_send_.begin() + i, 0);
    n_values_i_recv_offsets_[i] = std::accumulate(
        n_values_i_recv_.begin(), n_values_i_recv_.begin() + i, 0);
  }

  
  send_buffer_size_ =
      std::accumulate(n_values_i_send_.begin(), n_values_i_send_.end(), 0);
  recv_buffer_size_ =
      std::accumulate(n_values_i_recv_.begin(), n_values_i_recv_.end(), 0);
}

template <class idx_t>
idx_t Communicator<idx_t>::n_values_i_send(int mpi_rank) const {
  return n_values_i_send_[mpi_rank];
}

template <class idx_t>
idx_t Communicator<idx_t>::n_values_i_recv(int mpi_rank) const {
  return n_values_i_recv_[mpi_rank];
}

template <class idx_t> idx_t Communicator<idx_t>::send_buffer_size() const {
  return send_buffer_size_;
}
template <class idx_t> idx_t Communicator<idx_t>::recv_buffer_size() const {
  return recv_buffer_size_;
}

template class Communicator<uint32>;
template class Communicator<uint64>;

} // namespace mpi
} // namespace hydra
