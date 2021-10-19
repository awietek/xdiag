#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/mpi/alltoall.h>

namespace hydra::mpi {

class Communicator {
public:
  Communicator(std::vector<idx_t> const &n_values_i_send);

  std::vector<int> n_values_i_send() const { return n_values_i_send_; }
  std::vector<int> n_values_i_recv() const { return n_values_i_recv_; }

  std::vector<int> n_values_i_send_offsets() const {
    return n_values_i_send_offsets_;
  }
  std::vector<int> n_values_i_recv_offsets() const {
    return n_values_i_recv_offsets_;
  }

  idx_t send_buffer_size() const { return send_buffer_size_; }
  idx_t recv_buffer_size() const { return recv_buffer_size_; }

  template <class T>
  void add_to_send_buffer(int mpi_rank, T value, std::vector<T> &send_buffer) const {
    idx_t idx =
        n_values_i_send_offsets_[mpi_rank] + n_values_prepared_[mpi_rank];
    send_buffer[idx] = value;
    ++n_values_prepared_[mpi_rank];
  }

  void flush() {
    std::fill(n_values_prepared_.begin(), n_values_prepared_.end(), 0);
  }

  template <class T>
  void all_to_all(std::vector<T> const &send_buffer,
                  std::vector<T> &recv_buffer) const {
    Alltoallv<T>(send_buffer.data(), n_values_i_send_.data(),
                 n_values_i_send_offsets_.data(), recv_buffer.data(),
                 n_values_i_recv_.data(), n_values_i_recv_offsets_.data(),
                 MPI_COMM_WORLD);
  }

private:
  int mpi_rank_;
  int mpi_size_;

  mutable std::vector<int> n_values_prepared_;

  std::vector<int> n_values_i_send_;
  std::vector<int> n_values_i_recv_;
  std::vector<int> n_values_i_send_offsets_;
  std::vector<int> n_values_i_recv_offsets_;

  idx_t send_buffer_size_;
  idx_t recv_buffer_size_;
};

} // namespace hydra::mpi
