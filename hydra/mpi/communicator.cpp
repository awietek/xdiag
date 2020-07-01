#include "communicator.h"

#include <algorithm>

namespace hydra { namespace mpi {

Communicator::Communicator
(std::vector<int> const& n_values_i_send)
  : n_values_i_send_(n_values_i_send),
    n_values_i_recv_(n_values_i_send.size(), 0),
    n_values_i_send_offsets_(n_values_i_send.size(), 0),
    n_values_i_recv_offsets_(n_values_i_send.size(), 0),
    n_values_i_send_2_(n_values_i_send.size(), 0),
    n_values_i_recv_2_(n_values_i_send.size(), 0),
    n_values_i_send_offsets_2_(n_values_i_send.size(), 0),
    n_values_i_recv_offsets_2_(n_values_i_send.size(), 0)
{
  using std::accumulate;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
  assert(n_values_i_send.size() == mpi_size_);
  MPI_Alltoall(n_values_i_send.data(), 1, MPI_UNSIGNED_LONG, 
	       n_values_i_recv.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  for (int i=0; i<mpi_rank_; ++i)
    {
      n_values_i_send_offsets_ = accumulate(n_values_i_send_.begin(), 
					    n_values_i_send_.begin()+i, 0);
      n_values_i_recv_offsets_ = accumulate(n_values_i_recv_.begin(), 
					    n_values_i_recv_.begin()+i, 0);
      n_values_i_send_2_[i] = 2 * n_values_i_send_[i];
      n_values_i_recv_2_[i] = 2 * n_values_i_recv_[i];
      n_values_i_send_offsets_2_[i] = 2 * n_values_i_send_offsets_[i];
      n_values_i_recv_offsets_2_[i] = 2 * n_values_i_recv_offsets_[i];	
    }
  
  send_buffer_size_ = accumulate(n_values_i_send_.begin(), 
				 n_values_i_send_.end(), 0);
  recv_buffer_size_ = accumulate(n_values_i_recv_.begin(), 
				 n_values_i_recv_.end(), 0);
}

int Communicator::n_values_i_send(int mpi_rank) const
{ return n_values_i_send_[mpi_rank]; }
int Communicator::n_values_i_recv(int mpi_rank) const
{ return n_values_i_recv_[mpi_rank]; }
  
int Communicator::send_buffer_size() const
{ return send_buffer_size_; }
int Communicator::recv_buffer_size() const
{ return recv_buffer_size_; }

template <> void all_to_all(double* send_buffer, double* recv_buffer) const
{
  MPI_Alltoallv(send_buffer,
		n_values_i_send.data(),
		n_values_i_send_offsets_.data(),
		MPI_DOUBLE, recvbuf,
	        recv_buffer,
		n_values_i_recv.data(),
		n_values_i_recv_offsets_.data(),
		MPI_COMM_WORLD);
}

template <> void all_to_all(complex* send_buffer, complex* recv_buffer) const
{
  MPI_Alltoallv(send_buffer,
		n_values_i_send_2_.data(),
		n_values_i_send_offsets_2_.data(),
		MPI_DOUBLE, 
		n_values_i_recv_2_.data(),
		n_values_i_recv_offsets_2_.data(),
		MPI_DOUBLE,
		MPI_COMM_WORLD);
}
      

}  // namespace mpi
}  // namespace hydra
