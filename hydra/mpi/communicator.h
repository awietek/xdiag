#ifndef HYDRA_MPI_COMMUNCATOR_
#define HYDRA_MPI_COMMUNCATOR_

#include <mpi.h>
#include <vector>
#include <lila/allmpi.h>

#include <hydra/common.h>

namespace hydra {
namespace mpi {

template <class idx_t_> class Communicator {
public:
  using idx_t = idx_t_;
  
  Communicator(std::vector<idx_t> const &n_values_i_send);

  idx_t n_values_i_send(int mpi_rank) const;
  idx_t n_values_i_recv(int mpi_rank) const;

  idx_t send_buffer_size() const;
  idx_t recv_buffer_size() const;

  template <class T>
  void add_to_send_buffer(int mpi_rank, T value, T *send_buffer) {
    idx_t idx =
        n_values_i_send_offsets_[mpi_rank] + n_values_prepared_[mpi_rank];
    send_buffer[idx] = value;
    ++n_values_prepared_[mpi_rank];
  }

  void flush() {
    std::fill(n_values_prepared_.begin(), n_values_prepared_.end(), 0);
  }
  
  template <class T> void all_to_all(T *send_buffer, T *recv_buffer){
    lila::MPI_Alltoallv<T>(send_buffer, n_values_i_send_.data(),
                           n_values_i_send_offsets_.data(), recv_buffer,
                           n_values_i_recv_.data(),
                           n_values_i_recv_offsets_.data(), MPI_COMM_WORLD);
  }

private:
  int mpi_rank_;
  int mpi_size_;

  std::vector<int> n_values_prepared_;

  std::vector<int> n_values_i_send_;
  std::vector<int> n_values_i_recv_;
  std::vector<int> n_values_i_send_offsets_;
  std::vector<int> n_values_i_recv_offsets_;

  idx_t send_buffer_size_;
  idx_t recv_buffer_size_;
};

} // namespace mpi
} // namespace hydra

#endif
