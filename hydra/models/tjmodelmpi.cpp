#include "tjmodelmpi.h"
#include "hubbardmodeldetail.h"

namespace hydra { namespace models {

    template <class coeff_t, class state_t>
    TJModelMPI<coeff_t, state_t>::TJModelMPI
    (BondList bondlist, Couplings couplings, hilbertspaces::hubbard_qn qn)
      : n_sites_(bondlist.n_sites()),
	qn_(qn)
    {
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
     
      // Currently unused operatuors
      std::vector<std::pair<int, int>> currents_;
      std::vector<coeff_t> current_amplitudes_;
      std::vector<std::pair<int, int>> interactions_;
      std::vector<double> interaction_strengths_;
      double U_;

      // Use Hubbard routing so set interation terms
      hubbardmodeldetail::set_hubbard_terms<coeff_t>
      (bondlist, couplings, hoppings_, hopping_amplitudes_,
       currents_, current_amplitudes_,
       interactions_, interaction_strengths_,
       onsites_, onsite_potentials_,
       szszs_, szsz_amplitudes_,
       exchanges_, exchange_amplitudes_, U_);

      initialize();
    }
      
    template <class coeff_t, class state_t>
    void TJModelMPI<coeff_t, state_t>::apply_hamiltonian
    (const lila::VectorMPI<coeff_t>& in_vec, lila::VectorMPI<coeff_t>& out_vec,
     bool verbose)
    {
 
    }

    template <class coeff_t, class state_t>
    void TJModelMPI<coeff_t, state_t>::initialize()
    {
      hs_upspins_ = Spinhalf<state_t>(n_sites_, qn_.n_upspins);
      hs_downspins_ = Spinhalf<state_t>(n_sites_, qn_.n_downspins);
      hs_holes_ = Spinhalf<state_t>(n_sites_, 
				    n_sites_ - qn_.n_upspins);

      indexing_holes_= IndexTable<Spinhalf<state_t>, uint64>(hs_holes_);


      // clear previous upspin, downspin members
      my_upspins_.clear();
      my_upspins_offset_.clear();
      my_downspins_.clear();
      my_downspins_offset_.clear();

      // Collect all upspin configurations belonging to this mpi_rank
      uint64 offset = 0;
      for (auto state : hs_upspins_)
	{
	  if (mpi_rank_of_spins(state) == mpi_rank_) 
	    {
	      my_upspins_.push_back(state);
	      my_upspins_offset_[state] = offset;	
	      offset += hs_holes_.size();
	    }
	}

      // Check whether every upspin configuration belongs to some process
      uint64 dim_test = my_upspins_.size();
      MPI_Allreduce(MPI_IN_PLACE, &dim_test, 1, MPI_UNSIGNED_LONG, MPI_SUM, 
		    MPI_COMM_WORLD);
      assert(dim_test == hs_upspins_.size());

      // Set dimensions
      local_dim_ = offset;
      MPI_Allreduce(&local_dim_, &dim_, 1, MPI_UNSIGNED_LONG, MPI_SUM, 
		    MPI_COMM_WORLD);
      MPI_Allreduce(&local_dim_, &max_local_dim_, 1, MPI_UNSIGNED_LONG, MPI_MAX,
		    MPI_COMM_WORLD);
      MPI_Allreduce(&local_dim_, &min_local_dim_, 1, MPI_UNSIGNED_LONG, MPI_MIN,
		    MPI_COMM_WORLD);

      // Collect all downspin configurations belonging to this mpi_rank
      offset = 0;
      for (auto state : hs_downspins_)
	{
	  if (mpi_rank_of_spins(state) == mpi_rank_) 
	    {
	      my_downspins_.push_back(state);
	      my_downspins_offset_[state] = offset;		
	      offset += hs_holes_.size();
	    }
	}
      local_dim_downspins_ = offset;

      // Check whether also local downspin dimensions sum up to the total dim
      uint64 dim_downspins;
      MPI_Allreduce(&local_dim_downspins_, &dim_downspins, 1, MPI_UNSIGNED_LONG,
		    MPI_SUM, MPI_COMM_WORLD);
      assert(dim_downspins == dim_);


      ///////////////////////////////////////////////
      // Determine Communication patterns
      // for forward sending downspins
      n_downspins_i_send_forward_.resize(mpi_size_, 0);
      n_downspins_i_recv_forward_.resize(mpi_size_, 0);
      n_downspins_i_send_forward_offsets_.resize(mpi_size_, 0);
      n_downspins_i_recv_forward_offsets_.resize(mpi_size_, 0);

      for (state_t upspins : my_upspins_)
	for (state_t holes : hs_holes_)
	  {
	    state_t downspins = up_hole_to_down(upspins, holes);
	    int destination_mpi_rank = mpi_rank_of_spins(downspins);
	    ++n_downspins_i_send_forward_[destination_mpi_rank];
	  }
      
      // Communicate number of send/receive states
      MPI_Alltoall(n_downspins_i_send_forward_.data(), 1, MPI_INT,
		   n_downspins_i_recv_forward_.data(), 1, MPI_INT,
		   MPI_COMM_WORLD);

      // Compute offsets and total number of send/receive states
      for(int m = 0; m < mpi_size_; ++m)
	for(int n = 0; n < m; ++n)
	  {
	    n_downspins_i_send_forward_offsets_[m] += n_downspins_i_send_forward_[n];
	    n_downspins_i_recv_forward_offsets_[m] += n_downspins_i_recv_forward_[n];
	  }

      sum_n_downspins_i_send_forward_ =
	std::accumulate(n_downspins_i_send_forward_.begin(),
			n_downspins_i_send_forward_.end(), 0);
      sum_n_downspins_i_recv_forward_ =
	std::accumulate(n_downspins_i_recv_forward_.begin(),
			n_downspins_i_recv_forward_.end(), 0);

      // Communicate which downstates are sent/received by whom
      std::vector<state_t> downspins_i_send_forward_(sum_n_downspins_i_send_forward_, 0);
      std::vector<state_t> upspins_i_send_forward_(sum_n_downspins_i_send_forward_, 0);

      std::vector<uint64> n_downspins_already_prepared(mpi_size_, 0);
      for (state_t upspins : my_upspins_)
	for (state_t holes : hs_holes_)
	  {
	    state_t downspins = up_hole_to_down(upspins, holes);
	    int destination_mpi_rank = mpi_rank_of_spins(downspins);
	    
	    uint64 idx = n_downspins_i_send_forward_offsets_[destination_mpi_rank] +
		n_downspins_already_prepared[destination_mpi_rank]++;
	     downspins_i_send_forward_[idx] = downspins;
	     upspins_i_send_forward_[idx] = upspins;
	  }
      downspins_i_recv_forward_.resize(sum_n_downspins_i_recv_forward_);
      upspins_i_recv_forward_.resize(sum_n_downspins_i_recv_forward_);


      MPI_Alltoallv(downspins_i_send_forward_.data(), 
		    n_downspins_i_send_forward_.data(),
		    n_downspins_i_send_forward_offsets_.data(), 
		    MPI_UNSIGNED,
		    downspins_i_recv_forward_.data(), 		      
		    n_downspins_i_recv_forward_.data(),
		    n_downspins_i_recv_forward_offsets_.data(), 
		    MPI_UNSIGNED, MPI_COMM_WORLD);
      MPI_Alltoallv(upspins_i_send_forward_.data(), 
		    n_downspins_i_send_forward_.data(),
		    n_downspins_i_send_forward_offsets_.data(), 
		    MPI_UNSIGNED,
		    upspins_i_recv_forward_.data(), 		      
		    n_downspins_i_recv_forward_.data(),
		    n_downspins_i_recv_forward_offsets_.data(), 
		    MPI_UNSIGNED, MPI_COMM_WORLD);

      // downspins_i_send_forward_ and upspins_i_send_forward_ not neede anymore
      downspins_i_send_forward_.clear();
      upspins_i_send_forward_.clear();
      downspins_i_send_forward_.shrink_to_fit();
      upspins_i_send_forward_.shrink_to_fit();



      ///////////////////////////////////////////////
      // Determine Communication patterns
      // for back sending downspins
      n_upspins_i_send_back_.resize(mpi_size_, 0);
      n_upspins_i_recv_back_.resize(mpi_size_, 0);
      n_upspins_i_send_back_offsets_.resize(mpi_size_, 0);
      n_upspins_i_recv_back_offsets_.resize(mpi_size_, 0);


      for (state_t downspins : my_downspins_)
	for (state_t holes : hs_holes_)
	  {
	    state_t upspins = down_hole_to_up(downspins, holes);
	    int destination_mpi_rank = mpi_rank_of_spins(upspins);
	    ++n_downspins_i_send_back_[destination_mpi_rank];
	  }

      
      // Communicate number of send/receive states
      MPI_Alltoall(n_upspins_i_send_back_.data(), 1, MPI_INT,
		   n_upspins_i_recv_back_.data(), 1, MPI_INT, MPI_COMM_WORLD);

      
      // Compute offsets and total number of send/receive states
      for(int m = 0; m < mpi_size_; ++m)
	for(int n = 0; n < m; ++n)
	  {
	    n_upspins_i_send_back_offsets_[m] += n_upspins_i_send_back_[n];
	    n_upspins_i_recv_back_offsets_[m] += n_upspins_i_recv_back_[n];
	  }

      sum_n_upspins_i_send_back_ = std::accumulate(n_upspins_i_send_back_.begin(),
						   n_upspins_i_send_back_.end(), 0);
      sum_n_upspins_i_recv_back_ = std::accumulate(n_upspins_i_recv_back_.begin(),
						   n_upspins_i_recv_back_.end(), 0);

      
      // Communicate which upstates are sent/received by whom
      std::vector<state_t> downspins_i_send_back(sum_n_upspins_i_send_back_, 0);
      std::vector<state_t> upspins_i_send_back(sum_n_upspins_i_send_back_, 0);


      std::vector<uint64> n_upspins_already_prepared(mpi_size_, 0);
      for (state_t downspins : my_downspins_)
	for (state_t holes : hs_holes_)
	  {
	    state_t upspins = down_hole_to_up(downspins, holes);
	    int destination_mpi_rank = mpi_rank_of_spins(upspins);

	    uint64 idx = n_upspins_i_send_back_offsets_[destination_mpi_rank] +
		n_upspins_already_prepared[destination_mpi_rank]++; 
	      downspins_i_send_back[idx] = downspins;
	      upspins_i_send_back[idx] = upspins;
	  }
      
      downspins_i_recv_back_.resize(sum_n_upspins_i_recv_back_);
      upspins_i_recv_back_.resize(sum_n_upspins_i_recv_back_);
      
      MPI_Alltoallv(downspins_i_send_back_.data(), 
		    n_upspins_i_send_back_.data(),
		    n_upspins_i_send_back_offsets_.data(), 
		    MPI_UNSIGNED,
		    downspins_i_recv_back_.data(), 		      
		    n_upspins_i_recv_back_.data(),
		    n_upspins_i_recv_back_offsets_.data(), 
		    MPI_UNSIGNED, MPI_COMM_WORLD);
      MPI_Alltoallv(upspins_i_send_back_.data(), 
		    n_upspins_i_send_back_.data(),
		    n_upspins_i_send_back_offsets_.data(), 
		    MPI_UNSIGNED,
		    upspins_i_recv_back_.data(), 		      
		    n_upspins_i_recv_back_.data(),
		    n_upspins_i_recv_back_offsets_.data(), 
		    MPI_UNSIGNED, MPI_COMM_WORLD);

      
      // downspins_i_send_back_ and upspins_i_send_back_ not needed anymore
      downspins_i_send_back_.clear();
      upspins_i_send_back_.clear();
      downspins_i_send_back_.shrink_to_fit();
      upspins_i_send_back_.shrink_to_fit();
      

      // Allocate send/receive buffers
      buffer_size_ = std::max(std::max(sum_n_downspins_i_send_forward_, 
				       sum_n_downspins_i_recv_forward_),
			      std::max(sum_n_upspins_i_send_back_, 
				       sum_n_upspins_i_recv_back_));

      try 								       
	{
	  send_buffer_.resize(buffer_size_);
	  recv_buffer_.resize(buffer_size_);
	}
      catch(...)
	{
	  std::cerr << "[ " << mpi_rank_
		    << " ] Error: Could not allocate send/receive buffers!" 
		    << std::endl << std::flush;
	  MPI_Abort(MPI_COMM_WORLD,3);
	}
      
    }

    template class TJModelMPI<double, uint32>;
    template class TJModelMPI<complex, uint32>;
    template class TJModelMPI<double, uint64>;
    template class TJModelMPI<complex, uint64>;
       
  }  // namespace models
}  // namespace hydra
