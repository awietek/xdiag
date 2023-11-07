#include "communicator.h"
#include <mpi.h>
#include <vector>

#include <numeric>

namespace hydra::mpi {

Communicator::Communicator(std::vector<int64_t> const &n_values_i_send)
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

int64_t Communicator::n_values_i_send(int mpi_rank) const {
  return n_values_i_send_[mpi_rank];
}
int64_t Communicator::n_values_i_recv(int mpi_rank) const {
  return n_values_i_recv_[mpi_rank];
}

int64_t Communicator::n_values_i_send_offsets(int mpi_rank) const {
  return n_values_i_send_offsets_[mpi_rank];
}
int64_t Communicator::n_values_i_recv_offsets(int mpi_rank) const {
  return n_values_i_recv_offsets_[mpi_rank];
}

int64_t Communicator::send_buffer_size() const { return send_buffer_size_; }
int64_t Communicator::recv_buffer_size() const { return recv_buffer_size_; }

int64_t Communicator::n_values_prepared(int mpi_rank) const {
  return n_values_prepared_[mpi_rank];
}

template <class T>
void Communicator::add_to_send_buffer(int mpi_rank, T value,
                                      T *send_buffer) const {
  int64_t idx =
      n_values_i_send_offsets_[mpi_rank] + n_values_prepared_[mpi_rank];
  send_buffer[idx] = value;
  ++n_values_prepared_[mpi_rank];
}
template void
Communicator::add_to_send_buffer<double>(int, double value,
                                         double *send_buffer) const;
template void
Communicator::add_to_send_buffer<complex>(int, complex value,
                                          complex *send_buffer) const;

void Communicator::flush() {
  std::fill(n_values_prepared_.begin(), n_values_prepared_.end(), 0);
}

template <class T>
void Communicator::all_to_all(const T *send_buffer, T *recv_buffer) const {
  Alltoallv<T>(
      const_cast<T *>(send_buffer), const_cast<int *>(n_values_i_send_.data()),
      const_cast<int *>(n_values_i_send_offsets_.data()),
      const_cast<T *>(recv_buffer), const_cast<int *>(n_values_i_recv_.data()),
      const_cast<int *>(n_values_i_recv_offsets_.data()), MPI_COMM_WORLD);
}

template void Communicator::all_to_all<double>(const double *send_buffer,
                                               double *recv_buffer) const;
template void Communicator::all_to_all<complex>(const complex *send_buffer,
                                                complex *recv_buffer) const;

} // namespace hydra::mpi
