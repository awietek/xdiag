#include "tjmodelmpi.h"

#include <hydra/utils/bitops.h>
#include <hydra/utils/combinatorics.h>
#include "hubbardmodeldetail.h"

namespace hydra { namespace models {

    using combinatorics::down_hole_to_up;
    using combinatorics::up_hole_to_down;

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
      using utils::popcnt;
      using utils::gbit;
      using utils::gbits;
      using namespace hydra::combinatorics;

      verbose = true;

      Zeros(out_vec.vector_local());

      // Apply szsz terms
      
      int szsz_idx = 0;
      for (auto pair : szszs_)
      {
        const int s1 = std::min(pair.first, pair.second);
        const int s2 = std::max(pair.first, pair.second);
        const double jz = szsz_amplitudes_[szsz_idx]*0.25;
        if (std::abs(jz) > 1e-14)
          {
            uint64 upspin_idx = 0;
            for (const state_t& upspins : my_upspins_)
            {
            uint64 upspin_offset = my_upspins_offset_[upspins];
            uint64 downspin_offset=0;
            for (state_t holes: hs_holes_in_ups_)
              {
                state_t downspins = up_hole_to_down(upspins, holes);
                uint64 idx = upspin_offset + downspin_offset;
                auto coeff =
            jz*(((double)gbit(upspins, s1) - (double)gbit(downspins, s1)) *
          ((double)gbit(upspins, s2) - (double)gbit(downspins, s2)));
                out_vec(idx) += coeff * in_vec(idx);
                ++downspin_offset;
              }
            ++upspin_idx;
            }
          }
        ++szsz_idx;
      }
      


      // Spin exchange terms
      
      int exchange_idx=0;
      for (auto pair: exchanges_)
      {
        const int s1 = std::min(pair.first, pair.second);
        const int s2 = std::max(pair.first, pair.second);
        const coeff_t jx = exchange_amplitudes_[exchange_idx]*0.5;
        const state_t flipmask = ((state_t)1 << s1) | ((state_t)1 << s2);

        // Find out how many states is sent to each process
        std::vector<int> n_states_i_send(mpi_size_, 0);
        if (std::abs(jx) > 1e-14)
        {
          // Flip states and check out how much needs to be communicated
          for (const state_t& upspins : my_upspins_)
          for (state_t holes : hs_holes_in_ups_)
            {
              state_t downspins = up_hole_to_down(upspins, holes);

              if ((popcnt(upspins & flipmask) == 1) && 
                  (popcnt(downspins & flipmask) == 1))
                  {
                  state_t flipped_upspins = upspins ^ flipmask;
                  int target = mpi_rank_of_spins(flipped_upspins);
                  ++n_states_i_send[target];
                  }
            }

            // Exchange information on who sends how much to whom
            std::vector<int> n_states_i_recv(mpi_size_, 0);
            MPI_Alltoall(n_states_i_send.data(), 1, MPI_INT, n_states_i_recv.data(),
              1, MPI_INT, MPI_COMM_WORLD);

            // Sum up states sent/recvd and eventually resize buffers
            uint64 sum_n_states_i_send = std::accumulate(n_states_i_send.begin(), 
              n_states_i_send.end(), (uint64)0);
            uint64 sum_n_states_i_recv = std::accumulate(n_states_i_recv.begin(),
              n_states_i_recv.end(), (uint64)0);

            if (sum_n_states_i_send > send_buffer_.size())
              send_buffer_.resize(sum_n_states_i_send);
            if (sum_n_states_i_recv > recv_buffer_.size())
              recv_buffer_.resize(sum_n_states_i_recv);

            // Compute offsets of send/recieve states
            std::vector<int> n_states_i_send_offsets(mpi_size_, 0);
            std::vector<int> n_states_i_recv_offsets(mpi_size_, 0);
            for (int m = 0; m < mpi_size_; ++m) 
              for (int n = 0; n < m; ++n)
              {
                n_states_i_send_offsets[m] += n_states_i_send[n];
                n_states_i_recv_offsets[m] += n_states_i_recv[n];
              }
            
            // Flip states and check how much needs to be communicated
            std::vector<int> n_states_prepared(mpi_size_, 0);
            for (const state_t& upspins : my_upspins_)
            {
              uint64 upspin_offset = my_upspins_offset_[upspins];
              uint64 downspin_offset=0;
              for (state_t holes: hs_holes_in_ups_)
                {
                  state_t downspins = up_hole_to_down(upspins, holes);

                  if ((popcnt(upspins & flipmask) == 1) &&
                (popcnt(downspins & flipmask) == 1))
              {
                uint64 idx = upspin_offset + downspin_offset;
                state_t flipped_upspins = upspins ^ flipmask;
                int target = mpi_rank_of_spins(flipped_upspins);
                int send_idx = n_states_i_send_offsets[target] +
                  n_states_prepared[target];

                send_buffer_[send_idx] = in_vec.vector_local()(idx);
                ++n_states_prepared[target];
              }
                  ++downspin_offset;
              }
            }

                 //LilaPrint(in_vec.vector_local());
                 //for (int i=0; i<sum_n_states_i_send; ++i)
                 //	printf("send_buf[%d] = %f\n", i, send_buffer_[i]);


                // Check, whether correct number of states has been prepared
                for(int m = 0; m < mpi_size_; ++m)
            assert(n_states_prepared[m] == n_states_i_send[m]);

                // Alltoall call
                lila::MPI_Alltoallv<coeff_t>
            (send_buffer_.data(), n_states_i_send.data(),
             n_states_i_send_offsets.data(),
             recv_buffer_.data(), n_states_i_recv.data(),
             n_states_i_recv_offsets.data(),
             MPI_COMM_WORLD);

                // for (int i=0; i<sum_n_states_i_recv; ++i)
                // 	printf("recv_buf[%d] = %f\n", i, recv_buffer_[i]);

                // Get the original upspin configuration and its source proc
                std::vector<std::vector<state_t>> upspins_i_get_from_proc(mpi_size_);
                for (const state_t& upspins : my_upspins_)
            {
              state_t flipped_upspins = upspins ^ flipmask;
              int source = mpi_rank_of_spins(flipped_upspins);
              upspins_i_get_from_proc[source].push_back(upspins);
            }

                // for (int i=0; i<mpi_size_; ++i)
                // 	for (auto upspins : upspins_i_get_from_proc[i])
                // 	  printf("upspins_i_get: %d\n", upspins);


                // Sort the upspins and add coefficients to outvec
                //LilaPrint(out_vec.vector_local());

                uint64 recv_idx=0;
                for (int m=0; m<mpi_size_; ++m)
            {

              // Sort according to order of flipped upspins
              std::sort(upspins_i_get_from_proc[m].begin(),
                  upspins_i_get_from_proc[m].end(),
                  [&flipmask](state_t const& a, state_t const& b)
                  {
                    state_t flipped_a = a ^ flipmask;
                    state_t flipped_b = b ^ flipmask;
                    return flipped_a < flipped_b;
                  });

              for (const state_t& upspins : upspins_i_get_from_proc[m])
                {
                  uint64 upspin_offset = my_upspins_offset_[upspins];

                  uint64 downspin_offset=0;
                  for (state_t holes : hs_holes_in_ups_)
              {
                    state_t downspins = up_hole_to_down(upspins, holes);
                if ((popcnt(upspins & flipmask) == 1) &&
                    (popcnt(downspins & flipmask) == 1))
                  {
                    uint64 target_idx = upspin_offset + downspin_offset;
                     //printf("upspins_offset: %d, target_idx: %d, recv_idx: %d, coeff: %f\n",
                     	//     upspin_offset, target_idx, recv_idx, jx * recv_buffer_[recv_idx] );

                    out_vec.vector_local()(target_idx) += jx *
                recv_buffer_[recv_idx];
                // LilaPrint(out_vec.vector_local());
                    ++recv_idx;
                  }
                ++downspin_offset;
              }
              }
            }
                 //LilaPrint(out_vec.vector_local());

              }
            ++exchange_idx;
          }


      // Apply hoppings on downspins
      double t1 = MPI_Wtime();
      int hopping_idx=0;
      for (auto pair : hoppings_)
      	{
      	  const int s1 = std::min(pair.first, pair.second); 
      	  const int s2 = std::max(pair.first, pair.second);
      	  const coeff_t t = hopping_amplitudes_[hopping_idx];
      	  const uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);
      	  if (std::abs(t) > 1e-14)
      	    {
      	      // Loop over all configurations
      	      uint64 upspin_idx = 0;

      	      for (state_t upspins : my_upspins_)
      		{
      		  if ((upspins & flipmask) == 0)  // t-J hard-core constraint
      		    {
      		      uint64 upspin_offset = my_upspins_offset_[upspins];
      		      uint64 downspin_offset = 0;
      		      for (state_t holes : hs_holes_in_ups_)
      			{
      			  state_t downspins = up_hole_to_down(upspins, holes);

      			  // Check if hopping is possible
      			  if (((downspins & flipmask) != 0) && 
      			      ((downspins & flipmask) != flipmask))
      			    {
      			      double fermi = 
      				popcnt(gbits(downspins, s2-s1-1, s1+1)) % 2==0 ? 1. : -1.;

      			      uint64 idx = upspin_offset + downspin_offset;
      			      state_t new_downspins = downspins ^ flipmask;
			      state_t new_holes = up_down_to_hole(upspins, new_downspins);
      			      uint64 new_idx = upspin_offset + indexing_holes_in_ups_.index(new_holes);
      			      out_vec.vector_local()(new_idx) -= fermi * t * in_vec.vector_local()(idx);
      			    }

      			  ++downspin_offset;
      			}
      		    }
      		  ++upspin_idx;
      		}
      	    }
      	  ++hopping_idx;
      	}  // hopping on downspin
      double t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  down: %3.4f\n", t2-t1); 


      // Send configurations 
      t1 = MPI_Wtime();
      std::vector<uint64> n_states_already_prepared(mpi_size_, 0);
      uint64 downspin_offset=0;
      for (state_t holes : hs_holes_in_ups_)
      	{
      	  for (uint64 upspin_idx = 0; upspin_idx < my_upspins_.size(); ++upspin_idx)
      	    {
      	      state_t upspins = my_upspins_[upspin_idx];
      	      state_t downspins = up_hole_to_down(upspins, holes);
      	      int destination_mpi_rank = mpi_rank_of_spins(downspins);
      	      uint64 send_idx = n_downspins_i_send_forward_offsets_[destination_mpi_rank] +
      		n_states_already_prepared[destination_mpi_rank]++; 
      	      uint64 upspin_offset = my_upspins_offset_[upspins];
      	      uint64 idx = upspin_offset + downspin_offset;
      	      send_buffer_[send_idx] = in_vec.vector_local()(idx);
      	    }
      	  ++downspin_offset;
      	}
      t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  send forward (prepare): %3.4f\n", t2-t1); 


      t1 = MPI_Wtime();
      lila::MPI_Alltoallv<coeff_t>
      	(send_buffer_.data(), n_downspins_i_send_forward_.data(), 
      	 n_downspins_i_send_forward_offsets_.data(), 
      	 recv_buffer_.data(), n_downspins_i_recv_forward_.data(), 
      	 n_downspins_i_recv_forward_offsets_.data(), 
      	 MPI_COMM_WORLD);
      t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  send forward (Alltoallv): %3.4f\n", t2-t1); 


      // Sort the received data into the send buffer
      t1 = MPI_Wtime();
      std::fill(send_buffer_.begin(), send_buffer_.end(), 0);
      for (uint64 idx = 0; idx < sum_n_downspins_i_recv_forward_; ++idx)
      	{
      	  state_t upspins = upspins_i_recv_forward_[idx];
      	  state_t downspins = downspins_i_recv_forward_[idx];	    
	  state_t holes = down_up_to_hole(downspins, upspins);

      	  uint64 sorted_idx = my_downspins_offset_[downspins] + 
      	    indexing_holes_in_downs_.index(holes);
      	  send_buffer_[sorted_idx] = recv_buffer_[idx];
      	}
      t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  sort forward: %3.4f\n", t2-t1); 


      std::fill(recv_buffer_.begin(), recv_buffer_.end(), 0);


      // Apply hoppings on upspins
      t1 = MPI_Wtime();
      hopping_idx = 0;
      for (auto pair : hoppings_)
      	{
      	  const int s1 = std::min(pair.first, pair.second); 
      	  const int s2 = std::max(pair.first, pair.second);
      	  const coeff_t t = hopping_amplitudes_[hopping_idx];
      	  const uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);
      	  if (std::abs(t) > 1e-14)
      	    {
      	      // Loop over all configurations
      	      uint64 downspin_idx = 0;
      	      for (state_t downspins : my_downspins_)
      		{
      		  if ((downspins & flipmask) == 0)  // t-J hard-core constraint
      		    {
      		      uint64 downspin_offset = my_downspins_offset_[downspins];
      		      uint64 upspin_offset = 0;
      		      for (state_t holes : hs_holes_in_downs_)
      			{
      			  state_t upspins = down_hole_to_up(downspins, holes);

      			  // Check if hopping is possible
      			  if (((upspins & flipmask) != 0) && 
      			      ((upspins & flipmask) != flipmask))
      			    {
      			      double fermi = 
      				popcnt(gbits(upspins, s2-s1-1, s1+1)) % 2==0 ? 1. : -1.;
      			      uint64 idx = upspin_offset + downspin_offset;
      			      state_t new_upspins = upspins ^ flipmask;
      			      state_t new_holes = down_up_to_hole(downspins, new_upspins);
      			      uint64 new_idx = downspin_offset + indexing_holes_in_downs_.index(new_holes);
      			      recv_buffer_[new_idx] -= fermi * t * send_buffer_[idx];
      			    }

      			  ++upspin_offset;
      			}
      		    }
      		  ++downspin_idx;
      		}
      	    }
      	  ++hopping_idx;
      	}  // hopping on upspins
      t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  up:   %3.4f\n", t2-t1); 


      // Send back the resulting configurations 
      t1 = MPI_Wtime();
      std::fill(n_states_already_prepared.begin(), n_states_already_prepared.end(), 0);
      uint64 upspin_offset=0;
      for (state_t holes : hs_holes_in_downs_)
      	{
      	  for (uint64 downspin_idx = 0; downspin_idx < my_downspins_.size(); ++downspin_idx)
      	    {
      	      state_t downspins = my_downspins_[downspin_idx];
      	      state_t upspins = down_hole_to_up(downspins, holes);
      	      int destination_mpi_rank = mpi_rank_of_spins(upspins);
	      uint64 send_idx = n_upspins_i_send_back_offsets_[destination_mpi_rank] +
		n_states_already_prepared[destination_mpi_rank]++; 
      	      uint64 downspin_offset = my_downspins_offset_[downspins];
      	      uint64 idx = upspin_offset + downspin_offset;
      	      send_buffer_[send_idx] = recv_buffer_[idx];
      	    }
      	  ++upspin_offset;
      	}
      t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  send back (prepare):   %3.4f\n", t2-t1); 


      t1 = MPI_Wtime();
      lila::MPI_Alltoallv<coeff_t>
      	(send_buffer_.data(), n_upspins_i_send_back_.data(), 
      	 n_upspins_i_send_back_offsets_.data(), 
      	 recv_buffer_.data(), n_upspins_i_recv_back_.data(), 
      	 n_upspins_i_recv_back_offsets_.data(), 
      	 MPI_COMM_WORLD);
      t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  send back (Alltoall):   %3.4f\n", t2-t1); 


      // Sort the received data into the send buffer
      t1 = MPI_Wtime();
      std::fill(send_buffer_.begin(), send_buffer_.end(), 0);
      for (uint64 idx = 0; idx < sum_n_upspins_i_recv_back_; ++idx)
      	{
      	  state_t upspins = upspins_i_recv_back_[idx];
      	  state_t downspins = downspins_i_recv_back_[idx];
	  state_t holes = up_down_to_hole(upspins, downspins);
      	  uint64 sorted_idx = my_upspins_offset_[upspins] + 
      	    indexing_holes_in_downs_.index(holes);
      	  send_buffer_[sorted_idx] = recv_buffer_[idx];
      	}

      t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  sort back:   %3.4f\n", t2-t1); 


      // Fill in from send buffer
      t1 = MPI_Wtime();
      for (int k=0; k<in_vec.vector_local().size(); ++k)
      	out_vec.vector_local()(k) += send_buffer_[k];
      t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  fill outvec:   %3.4f\n", t2-t1); 

    }

    template <class coeff_t, class state_t>
    void TJModelMPI<coeff_t, state_t>::initialize()
    {
      hs_upspins_ = Spinhalf<state_t>(n_sites_, qn_.n_upspins);
      hs_downspins_ = Spinhalf<state_t>(n_sites_, qn_.n_downspins);
      hs_holes_in_ups_ = Spinhalf<state_t>(n_sites_ - qn_.n_upspins, 
					   qn_.n_downspins);
      hs_holes_in_downs_ = Spinhalf<state_t>(n_sites_ - qn_.n_downspins, 
					     qn_.n_upspins);

      indexing_holes_in_ups_= LinTable<Spinhalf<state_t>, uint64>(hs_holes_in_ups_);
      indexing_holes_in_downs_= LinTable<Spinhalf<state_t>, uint64>(hs_holes_in_downs_);

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
	      offset += hs_holes_in_ups_.size();
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
	      offset += hs_holes_in_downs_.size();
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
	{
	  // if (mpi_rank_ == 0)
	  //   std::cout << "up " << PrintSpinhalf(n_sites_, upspins) << std::endl;
	
	  for (state_t holes : hs_holes_in_ups_)
	    {
	      state_t downspins = up_hole_to_down(upspins, holes);
	      // if (mpi_rank_ == 0)
	      // 	{
	      // 	  std::cout << "hl " << PrintSpinhalf(n_sites_, holes) << std::endl;
	      // 	  std::cout << "dn " << PrintSpinhalf(n_sites_, downspins) << std::endl;
	      // 	}

	      int destination_mpi_rank = mpi_rank_of_spins(downspins);
	      ++n_downspins_i_send_forward_[destination_mpi_rank];
	    }
	  // if (mpi_rank_ == 0) std::cout << std::endl;
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
	for (state_t holes : hs_holes_in_ups_)
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

      // downspins_i_send_forward_ and upspins_i_send_forward_ not needed anymore
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
	for (state_t holes : hs_holes_in_downs_)
	  {
	    state_t upspins = down_hole_to_up(downspins, holes);
	    int destination_mpi_rank = mpi_rank_of_spins(upspins);
	    ++n_upspins_i_send_back_[destination_mpi_rank];
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
      std::vector<state_t> downspins_i_send_back_(sum_n_upspins_i_send_back_, 0);
      std::vector<state_t> upspins_i_send_back_(sum_n_upspins_i_send_back_, 0);


      std::vector<uint64> n_upspins_already_prepared(mpi_size_, 0);
	  // if (mpi_rank_ == 0) std::cout << std::endl;

      for (state_t downspins : my_downspins_)
	{
	  // if (mpi_rank_ == 0)
	  //   std::cout << "dn " << PrintSpinhalf(n_sites_, downspins) << std::endl;
	
	for (state_t holes : hs_holes_in_downs_)
	  {
	    state_t upspins = down_hole_to_up(downspins, holes);

	    // if (mpi_rank_ == 0)
	    //   {
	    // 	std::cout << "hl " << PrintSpinhalf(n_sites_, holes) << std::endl;
	    // 	std::cout << "up " << PrintSpinhalf(n_sites_, upspins) << std::endl;
	    //   }


	    int destination_mpi_rank = mpi_rank_of_spins(upspins);

	    uint64 idx = n_upspins_i_send_back_offsets_[destination_mpi_rank] +
		n_upspins_already_prepared[destination_mpi_rank]++; 
	      downspins_i_send_back_[idx] = downspins;
	      upspins_i_send_back_[idx] = upspins;
	  }
	  // if (mpi_rank_ == 0) std::cout << std::endl;

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
