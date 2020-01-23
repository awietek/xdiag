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
				    n_sites_ - qn_.n_upspins - qn_.n_downspins);

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

    }

    template class TJModelMPI<double, uint32>;
    template class TJModelMPI<complex, uint32>;
    template class TJModelMPI<double, uint64>;
    template class TJModelMPI<complex, uint64>;
       
  }  // namespace models
}  // namespace hydra
