#ifndef HYDRA_MPI_COMMUNCATOR_
#define HYDRA_MPI_COMMUNCATOR_

#include <mpi.h>
#include <vector>

namespace hydra { namespace mpi {

    class Communicator
    {
      Communicator(std::vector<uint64> const& n_values_i_send);

      uint64 n_values_i_send(int mpi_rank) const;
      uint64 n_values_i_recv(int mpi_rank) const;

      uint64 send_buffer_size() const;
      uint64 recv_buffer_size() const;

      template <class T>
      void all_to_all(T* send_buffer, T* recv_buffer) const;
      
    private:
      int mpi_rank_;
      int mpi_size_;
      
      std::vector<uint64> n_values_i_send_(mpi_size_, 0);
      std::vector<uint64> n_values_i_recv_(mpi_size_, 0);
      std::vector<uint64> n_values_i_send_offsets_(mpi_size_, 0);
      std::vector<uint64> n_values_i_recv_offsets_(mpi_size_, 0);

      // For complex values, same as above but multiplied with 2
      std::vector<uint64> n_values_i_send_2_(mpi_size_, 0);
      std::vector<uint64> n_values_i_recv_2_(mpi_size_, 0);
      std::vector<uint64> n_values_i_send_offsets_2_(mpi_size_, 0);
      std::vector<uint64> n_values_i_recv_offsets_2_(mpi_size_, 0);
      
      uint64 send_buffer_size_;
      uint64 recv_buffer_size_;
    }

  }  // namespace mpi
}  // namespace hydra

#endif
