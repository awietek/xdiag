#pragma once
#ifdef HYDRA_USE_MPI

#include <hydra/common.h>
#include <hydra/parallel/mpi/alltoall.h>
#include <hydra/parallel/mpi/buffer.h>
#include <mpi.h>

namespace hydra::mpi {

class Communicator {
public:
  Communicator() = default;
  Communicator(std::vector<int64_t> const &n_values_i_send);

  int64_t n_values_i_send(int mpi_rank) const;
  int64_t n_values_i_recv(int mpi_rank) const;

  int64_t n_values_i_send_offsets(int mpi_rank) const;
  int64_t n_values_i_recv_offsets(int mpi_rank) const;

  int64_t send_buffer_size() const;
  int64_t recv_buffer_size() const;

  int64_t n_values_prepared(int mpi_rank) const;

  template <class T>
  void add_to_send_buffer(int mpi_rank, T value, T *send_buffer) const;

  void flush();

  template <class T>
  void all_to_all(const T *send_buffer, T *recv_buffer) const;

private:
  int mpi_rank_;
  int mpi_size_;

  mutable std::vector<int> n_values_prepared_;

  std::vector<int> n_values_i_send_;
  std::vector<int> n_values_i_recv_;
  std::vector<int> n_values_i_send_offsets_;
  std::vector<int> n_values_i_recv_offsets_;

  int64_t send_buffer_size_;
  int64_t recv_buffer_size_;
};

} // namespace hydra::mpi
#endif
