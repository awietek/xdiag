#pragma once
#ifdef HYDRA_ENABLE_MPI

#include <mpi.h>
#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/mpi/alltoall.h>
#include <hydra/mpi/buffer.h>

namespace hydra::mpi {

class Communicator {
public:
  Communicator(std::vector<int64_t> const &n_values_i_send);

  std::vector<int> n_values_i_send() const { return n_values_i_send_; }
  std::vector<int> n_values_i_recv() const { return n_values_i_recv_; }

  std::vector<int> n_values_i_send_offsets() const {
    return n_values_i_send_offsets_;
  }
  std::vector<int> n_values_i_recv_offsets() const {
    return n_values_i_recv_offsets_;
  }

  int64_t send_buffer_size() const { return send_buffer_size_; }
  int64_t recv_buffer_size() const { return recv_buffer_size_; }

  template <class T>
  void add_to_send_buffer(int mpi_rank, T value, T *send_buffer) const {
    int64_t idx =
        n_values_i_send_offsets_[mpi_rank] + n_values_prepared_[mpi_rank];

    // if (typeid(T) == typeid(double)){
    // lila::Log("double idx {}", idx);
    // } else {
    //       lila::Log("cplx idx {}", idx);
    // }

    send_buffer[idx] = value;
    ++n_values_prepared_[mpi_rank];
  }

  void flush() {
    std::fill(n_values_prepared_.begin(), n_values_prepared_.end(), 0);
  }

  template <class T>
  void all_to_all(const T *send_buffer, T *recv_buffer) const {
    Alltoallv<T>(const_cast<T*>(send_buffer), const_cast<int*>(n_values_i_send_.data()),
                 const_cast<int*>(n_values_i_send_offsets_.data()), buffer.recv<T>(),
                 const_cast<int*>(n_values_i_recv_.data()),
		 const_cast<int*>(n_values_i_recv_offsets_.data()),
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

} // namespace hydra::mpi
#endif
