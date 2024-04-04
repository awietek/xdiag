#pragma once
#ifdef XDIAG_USE_MPI

#include <mpi.h>

#include <xdiag/common.hpp>
#include <xdiag/parallel/mpi/alltoall.hpp>
#include <xdiag/parallel/mpi/buffer.hpp>

namespace xdiag::mpi {

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

  void flush();

  template <class T>
  inline void add_to_send_buffer(int mpi_rank, T value, T *send_buffer) const {
    int64_t idx =
        n_values_i_send_offsets_[mpi_rank] + n_values_prepared_[mpi_rank];
    send_buffer[idx] = value;
    ++n_values_prepared_[mpi_rank];
  }

  template <class T>
  inline void all_to_all(const T *send_buffer, T *recv_buffer) const {
    Alltoallv<T>(const_cast<T *>(send_buffer),
                 const_cast<int *>(n_values_i_send_.data()),
                 const_cast<int *>(n_values_i_send_offsets_.data()),
                 const_cast<T *>(recv_buffer),
                 const_cast<int *>(n_values_i_recv_.data()),
                 const_cast<int *>(n_values_i_recv_offsets_.data()),
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

  int64_t send_buffer_size_;
  int64_t recv_buffer_size_;
};

} // namespace xdiag::mpi
#endif
