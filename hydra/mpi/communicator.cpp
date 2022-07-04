#include "communicator.h"
#include <mpi.h>
#include <vector>

#include <numeric>

namespace hydra::mpi {

Communicator::Communicator(std::vector<idx_t> const &n_values_i_send)
    : n_values_prepared_(n_values_i_send.size(), 0),
      n_values_i_recv_(n_values_i_send.size(), 0),
      n_values_i_send_offsets_(n_values_i_send.size(), 0),
      n_values_i_recv_offsets_(n_values_i_send.size(), 0) {
  for (auto num : n_values_i_send)
    n_values_i_send_.push_back((int)num);

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
  assert((int)n_values_i_send.size() == mpi_size_);
  MPI_Alltoall(n_values_i_send_.data(), 1, MPI_INT, n_values_i_recv_.data(), 1,
               MPI_INT, MPI_COMM_WORLD);

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

} // namespace hydra::mpi
