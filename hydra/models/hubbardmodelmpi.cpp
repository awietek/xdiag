#include <iostream>
#include <hydra/utils/bitops.h>
#include <hydra/utils/complex.h>
#include "hubbardmodeldetail.h"
#include "hubbardmodelmpi.h"

namespace hydra { namespace models {
    
    using hilbertspaces::Spinhalf;
    using hilbertspaces::PrintSpinhalf;
    using hilbertspaces::hubbard_qn;
    using indexing::IndexTable;

    template <class coeff_t, class state_t>
    HubbardModelMPI<coeff_t, state_t>::HubbardModelMPI
    (BondList bondlist, Couplings couplings, hilbertspaces::hubbard_qn qn)
      : n_sites_(bondlist.n_sites()),
	qn_(qn)
    {
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
     
      hubbardmodeldetail::set_hubbard_terms<coeff_t>
      (bondlist, couplings, hoppings_, hopping_amplitudes_,
       currents_, current_amplitudes_, interactions_, interaction_strengths_,
       onsites_, onsite_potentials_, U_);

      initialize();
    }

    template <class coeff_t, class state_t>
    void HubbardModelMPI<coeff_t, state_t>::set_qn(hilbertspaces::hubbard_qn qn)
    { 
      qn_ = qn; 
      initialize();
    }
      
    template <class coeff_t, class state_t>
    void HubbardModelMPI<coeff_t, state_t>::apply_hamiltonian
    (const lila::VectorMPI<coeff_t>& in_vec, lila::VectorMPI<coeff_t>& out_vec,
     bool verbose)
    {
      using utils::popcnt;
      using utils::gbit;
      using utils::gbits;

      if ((mpi_rank_== 0) && verbose)
	{
	  printf("max_local_dim: %lu\n", max_local_dim_);
	  printf("min_local_dim: %lu\n", min_local_dim_);
	}
      auto hs_downspins = Spinhalf<state_t>(n_sites_, qn_.n_downspins);
      Zeros(out_vec.vector_local());

      double t1 = MPI_Wtime();

      // Apply Diagonal U terms
      if (std::abs(U_) > 1e-14)
	{
	  uint64 upspin_idx = 0;
	  for (const state_t& upspins : my_upspins_)  // loop over upspins of process
	    {
	      uint64 upspin_offset = my_upspins_offset_[upspins];
	      uint64 downspin_offset=0;
	      for (state_t downspins : hs_downspins)  // loop over all downspins
		{
		  uint64 idx = upspin_offset + downspin_offset;
		  double coeff = U_*(double)popcnt(upspins & downspins);
		  out_vec.vector_local()(idx) += coeff*in_vec.vector_local()(idx);
		  ++downspin_offset;
		}
	      ++upspin_idx;
	    }
	}

     // Apply Diagonal V terms
      int interaction_idx=0;
      for (auto pair : interactions_)
	{
	  const int s1 = pair.first; 
	  const int s2 = pair.second;
	  const double V = interaction_strengths_[interaction_idx];

	  if (std::abs(V) > 1e-14)
	    {
	      uint64 upspin_idx = 0;
	      for (const state_t& upspins : my_upspins_)  // loop over upspins of process
		{
		  uint64 upspin_offset = my_upspins_offset_[upspins];
		  uint64 downspin_offset=0;
		  for (state_t downspins : hs_downspins)  // loop over all downspins
		    {
		      uint64 idx = upspin_offset + downspin_offset;
		      auto coeff = 	
			V * (double)((gbit(upspins, s1) + gbit(downspins, s1))*
				     (gbit(upspins, s2) + gbit(downspins, s2)));
		      out_vec(idx) += coeff * in_vec(idx); 
		      ++downspin_offset;
		    }
		  ++upspin_idx;
		}
	    }  // (std::abs(V) > 1e-14)

	  ++interaction_idx;
	}


     // Apply onsite chemical potential
      int onsite_idx=0;
      for (auto site : onsites_)
	{
	  const double mu = onsite_potentials_[onsite_idx];

	  if (std::abs(mu) > 1e-14)
	    {
	      uint64 upspin_idx = 0;
	      for (const state_t& upspins : my_upspins_)  // loop over upspins of process
		{
		  uint64 upspin_offset = my_upspins_offset_[upspins];
		  uint64 downspin_offset=0;
		  for (state_t downspins : hs_downspins)  // loop over all downspins
		    {
		      uint64 idx = upspin_offset + downspin_offset;
		      auto coeff = 	
			mu * (double)((gbit(upspins, site) + gbit(downspins, site)));
		      out_vec(idx) -= coeff * in_vec(idx); 
		      ++downspin_offset;
		    }
		  ++upspin_idx;
		}
	    }  // (std::abs(mu) > 1e-14)

	  ++onsite_idx;
	}


      double t2 = MPI_Wtime();
      if ((mpi_rank_== 0) && verbose) printf("  diag: %3.4f\n", t2-t1); 


      
      // Apply hoppings on downspins
      t1 = MPI_Wtime();
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
	      for (const state_t& upspins : my_upspins_) 
		{
		  uint64 upspin_offset = my_upspins_offset_[upspins];
		  uint64 downspin_offset = 0;
		  for (state_t downspins : hs_downspins)
		    {

		      // Check if hopping is possible
		      if (((downspins & flipmask) != 0) && 
			  ((downspins & flipmask) != flipmask))
			{
			  double fermi = 
			    popcnt(gbits(downspins, s2-s1-1, s1+1)) % 2==0 ? 1. : -1.;

			  uint64 idx = upspin_offset + downspin_offset;
			  state_t new_downspins = downspins ^ flipmask;
			  uint64 new_idx = upspin_offset + indexing_downspins_.index(new_downspins);
			  out_vec.vector_local()(new_idx) -= fermi * t * in_vec.vector_local()(idx);
			}

		      ++downspin_offset;
		    }
		  ++upspin_idx;
		}
	    }
	  ++hopping_idx;
	}  // hopping on downspin

      // Apply currents on downspins
      t1 = MPI_Wtime();
      int current_idx=0;
      for (auto pair : currents_)
	{
	  const int s1 = std::min(pair.first, pair.second); 
	  const int s2 = std::max(pair.first, pair.second);
	  const coeff_t t = current_amplitudes_[current_idx];
	  const uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);
	  // printf("down p1: %d, p2: %d, s1: %d, s2: %d\n",
	  // 	 pair.first, pair.second, s1, s2);
			  
	  if (std::abs(t) > 1e-14)
	    {
	      // Loop over all configurations
	      uint64 upspin_idx = 0;
	      for (const state_t& upspins : my_upspins_) 
		{
		  uint64 upspin_offset = my_upspins_offset_[upspins];
		  uint64 downspin_offset = 0;
		  for (state_t downspins : hs_downspins)
		    {

		      // Check if current is possible
		      if (((downspins & flipmask) != 0) && 
			  ((downspins & flipmask) != flipmask))
			{
			  double fermi =
			    popcnt(gbits(downspins, s2-s1-1, s1+1)) % 2==0 ? 1. : -1.;
			  double dir = gbit(downspins, pair.first) ? 1. : -1.;  // use pair, not s1 s2
			  // printf("down p1: %d, p2: %d, s1: %d, s2: %d, v1: %d, v2: %d, dir: %f\n",
			  // 	 pair.first, pair.second, s1, s2, gbit(downspins, pair.first), gbit(downspins, pair.second), dir);
			  
			  uint64 idx = upspin_offset + downspin_offset;
			  state_t new_downspins = downspins ^ flipmask;
			  uint64 new_idx = upspin_offset + indexing_downspins_.index(new_downspins);
			  out_vec.vector_local()(new_idx) -= fermi * dir * t * in_vec.vector_local()(idx);
			}

		      ++downspin_offset;
		    }
		  ++upspin_idx;
		}
	    }
	  ++hopping_idx;
	}  // hopping on downspin

      t2 = MPI_Wtime();
      // if ((mpi_rank_ == 0) && verbose) printf("  down: %3.4f\n", t2-t1); 


      //
      // Communication for hoppings on upspins
      //
	
      // Send configurations 
      t1 = MPI_Wtime();
      std::vector<uint64> n_states_already_prepared(mpi_size_, 0);
      uint64 downspin_offset=0;
      for (state_t downspins : hs_downspins)
	{	    
	  int destination_mpi_rank = mpi_rank_of_spins(downspins);
	  for (uint64 upspin_idx = 0; upspin_idx < my_upspins_.size(); ++upspin_idx)
	    {
	      uint64 send_idx = n_downspins_i_send_forward_offsets_[destination_mpi_rank] +
		n_states_already_prepared[destination_mpi_rank]++; 
	      uint64 upspin_offset = my_upspins_offset_[my_upspins_[upspin_idx]];
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


      // // DEBUG START
      // // Fat DEBUG Print
      // MPI_Barrier(MPI_COMM_WORLD);
	
      // for(int nt = 0; nt < mpi_size_; ++nt) 
      //   {
      //     if(mpi_rank_ == nt) 
      //       {
      // 	printf("[%d]:\n", mpi_rank_);

      // 	uint64 upspin_idx = 0;
      // 	for (const state_t& upspins : my_upspins_)  // loop over upspins of process
      // 	  {
      // 	    uint64 upspin_offset = my_upspins_offset_[upspins];
      // 	    uint64 downspin_offset=0;
      // 	    for (state_t downspins : hs_downspins)  // loop over all downspins
      // 	      {
      // 		uint64 idx = upspin_offset + downspin_offset;
      // 		printf("%s; %s %f -> [%d]\n", 
      // 		       PrintSpinhalf(n_sites_, upspins).c_str(),
      // 		       PrintSpinhalf(n_sites_, downspins).c_str(),
      // 		       in_vec.vector_local()(idx),
      // 		       mpi_rank_of_spins(downspins));
      // 		++downspin_offset;
      // 	      }
      // 	    ++upspin_idx;
      // 	  }

      // 	printf("\n");
      //       }
      //     MPI_Barrier(MPI_COMM_WORLD);
      //   }
      // MPI_Barrier(MPI_COMM_WORLD);
      // if(mpi_rank_ == 0) printf("SEND BUFFER --------------------------------\n");
      // MPI_Barrier(MPI_COMM_WORLD);
      // for(int nt = 0; nt < mpi_size_; ++nt) 
      //   {
      //     if(mpi_rank_ == nt) 
      //       {
      // 	printf("[%d]:\n", mpi_rank_);

      // 	for (int sendto = 0; sendto < mpi_size_; ++sendto)
      // 	  {
      // 	    int from = n_downspins_i_send_forward_offsets_[sendto];
      // 	    int to = (sendto == mpi_size_-1) ? sum_n_downspins_i_send_forward_ : n_downspins_i_send_forward_offsets_[sendto+1];
      // 	    for (int idx = from; idx < to; ++idx)
      // 	      printf("%f -> [%d] %s; %s\n", 
      // 		     send_buffer_[idx], sendto,
      // 		     PrintSpinhalf(n_sites_, upspins_i_send_forward_[idx]).c_str(),
      // 		     PrintSpinhalf(n_sites_, downspins_i_send_forward_[idx]).c_str());
      // 	  }
      // 	printf("\n");
      //       }
      //     MPI_Barrier(MPI_COMM_WORLD);
      //   }
      // MPI_Barrier(MPI_COMM_WORLD);



      // if(mpi_rank_ == 0) printf("RECV BUFFER --------------------------------\n");
      // MPI_Barrier(MPI_COMM_WORLD);
      // for(int nt = 0; nt < mpi_size_; ++nt) 
      //   {
      //     if(mpi_rank_ == nt) 
      //       {
      // 	printf("[%d]:\n", mpi_rank_);
      // 	for (int recvfrom = 0; recvfrom < mpi_size_; ++recvfrom)
      // 	  {
      // 	    int from = n_downspins_i_recv_forward_offsets_[recvfrom];
      // 	    int to = (recvfrom == mpi_size_-1) ? sum_n_downspins_i_recv_forward_ : n_downspins_i_recv_forward_offsets_[recvfrom+1];
      // 	    for (int idx = from; idx < to; ++idx)
      // 	      printf("%f <- [%d]  %s; %s\n", recv_buffer_[idx], recvfrom,
      // 		     PrintSpinhalf(n_sites_, upspins_i_recv_forward_[idx]).c_str(),
      // 		     PrintSpinhalf(n_sites_, downspins_i_recv_forward_[idx]).c_str());
		    
      // 	  }
      // 	printf("\n");
      //       }
      //     MPI_Barrier(MPI_COMM_WORLD);
      //   }
      // MPI_Barrier(MPI_COMM_WORLD);
      // // DEBUG END



      // Sort the received data into the send buffer
      t1 = MPI_Wtime();
      std::fill(send_buffer_.begin(), send_buffer_.end(), 0);
      for (uint64 idx = 0; idx < sum_n_downspins_i_recv_forward_; ++idx)
	{
	  state_t upspins = upspins_i_recv_forward_[idx];
	  state_t downspins = downspins_i_recv_forward_[idx];	    
	    
	  uint64 sorted_idx = my_downspins_offset_[downspins] + 
	    indexing_upspins_.index(upspins);
	  send_buffer_[sorted_idx] = recv_buffer_[idx];
	}
      t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  sort forward: %3.4f\n", t2-t1); 



      // // DEBUG START
      // if(mpi_rank_ == 0) printf("SORTED SEND BUFFER --------------------------------\n");
      // MPI_Barrier(MPI_COMM_WORLD);
      // for(int nt = 0; nt < mpi_size_; ++nt) 
      //   {
      //     if(mpi_rank_ == nt) 
      //       {
      // 	printf("[%d]:\n", mpi_rank_);
      // 	uint64 idx = 0;
      // 	for (const state_t& downspins : my_downspins_) 
      // 	  {
      // 	    for (state_t upspins : hs_upspins_)
      // 	      {
      // 		printf("%s; %s  %f\n", 
      // 		       PrintSpinhalf(n_sites_, upspins).c_str(),
      // 		       PrintSpinhalf(n_sites_, downspins).c_str(),
      // 		       send_buffer_[idx]);
      // 		++idx;
      // 	      }
      // 	  }
      // 	printf("\n");
      //       }
      //     MPI_Barrier(MPI_COMM_WORLD);
      //   }
      // MPI_Barrier(MPI_COMM_WORLD);
      // // DEBUG END


	
      // Apply upspin hoppings on the sorted send buffer
      // write the results to the receive buffer
	
      // // DEBUG START
      // // JUST COPY INVEC FOR DEBUGGING
      // std::copy(send_buffer_.begin(), send_buffer_.end(), recv_buffer_.begin());
      // // DEBUG END

      std::fill(recv_buffer_.begin(), recv_buffer_.end(), 0);

      t1 = MPI_Wtime();

      // Apply hoppings on upspins
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
	      for (const state_t& downspins : my_downspins_) 
		{
		  uint64 downspin_offset = my_downspins_offset_[downspins];
		  uint64 upspin_offset = 0;
		  for (state_t upspins : hs_upspins_)
		    {
		      // Check if hopping is possible
		      if (((upspins & flipmask) != 0) && 
			  ((upspins & flipmask) != flipmask))
			{
			  double fermi = 
			    popcnt(gbits(upspins, s2-s1-1, s1+1)) % 2==0 ? 1. : -1.;
			  uint64 idx = upspin_offset + downspin_offset;
			  state_t new_upspins = upspins ^ flipmask;
			  uint64 new_idx = downspin_offset + indexing_upspins_.index(new_upspins);
			  recv_buffer_[new_idx] -= fermi * t * send_buffer_[idx];
			}

		      ++upspin_offset;
		    }
		  ++downspin_idx;
		}
	    }
	  ++hopping_idx;
	}  // hopping on upspins


      // Apply currents on upspins
      current_idx = 0;
      for (auto pair : currents_)
	{
	  const int s1 = std::min(pair.first, pair.second); 
	  const int s2 = std::max(pair.first, pair.second);
	  const coeff_t t = current_amplitudes_[current_idx];
	  const uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);
	  // printf("up p1: %d, p2: %d, s1: %d, s2: %d\n",
	  // 	 pair.first, pair.second, s1, s2);
	  if (std::abs(t) > 1e-14)
	    {
	      // Loop over all configurations
	      uint64 downspin_idx = 0;
	      for (const state_t& downspins : my_downspins_) 
		{
		  uint64 downspin_offset = my_downspins_offset_[downspins];
		  uint64 upspin_offset = 0;
		  for (state_t upspins : hs_upspins_)
		    {
		      // Check if current is possible
		      if (((upspins & flipmask) != 0) && 
			  ((upspins & flipmask) != flipmask))
			{
			  double fermi = popcnt(gbits(upspins, s2-s1-1, s1+1)) % 2==0 ? 1. : -1.;
			  double dir = gbit(upspins, pair.first) ? 1. : -1.;  // use pair, not s1 s2
			  // printf("up p1: %d, p2: %d, s1: %d, s2: %d, v1: %d, v2: %d, dir: %f\n",
			  // 	 pair.first, pair.second, s1, s2, gbit(upspins, pair.first), gbit(upspins, pair.second), dir);
			  
			  uint64 idx = upspin_offset + downspin_offset;
			  state_t new_upspins = upspins ^ flipmask;
			  uint64 new_idx = downspin_offset + indexing_upspins_.index(new_upspins);
			  recv_buffer_[new_idx] -= fermi * dir * t * send_buffer_[idx];
			}

		      ++upspin_offset;
		    }
		  ++downspin_idx;
		}
	    }
	  ++current_idx;
	}  // current on upspins

      t2 = MPI_Wtime();
      // if ((mpi_rank_ == 0) && verbose) printf("  up:   %3.4f\n", t2-t1); 
	


      // // DEBUG START
      // if(mpi_rank_ == 0) printf("RESULTING RECV BUFFER --------------------------------\n");
      // MPI_Barrier(MPI_COMM_WORLD);
      // for(int nt = 0; nt < mpi_size_; ++nt) 
      //   {
      //     if(mpi_rank_ == nt) 
      //       {
      // 	printf("[%d]:\n", mpi_rank_);
      // 	uint64 idx = 0;
      // 	for (const state_t& downspins : my_downspins_) 
      // 	  {
      // 	    for (state_t upspins : hs_upspins_)
      // 	      {
      // 		printf("%s; %s  %f -> [%d]\n", 
      // 		       PrintSpinhalf(n_sites_, upspins).c_str(),
      // 		       PrintSpinhalf(n_sites_, downspins).c_str(),
      // 		       recv_buffer_[idx], mpi_rank_of_spins(upspins));
      // 		++idx;
      // 	      }
      // 	  }
      // 	printf("\n");
      //       }
      //     MPI_Barrier(MPI_COMM_WORLD);
      //   }
      // MPI_Barrier(MPI_COMM_WORLD);
      // // DEBUG END


	
      // Send back the resulting configurations 
      t1 = MPI_Wtime();
      std::fill(n_states_already_prepared.begin(), n_states_already_prepared.end(), 0);
      uint64 upspin_offset=0;
      for (state_t upspins : hs_upspins_)
	{	    
	  int destination_mpi_rank = mpi_rank_of_spins(upspins);
	  for (uint64 downspin_idx = 0; downspin_idx < my_downspins_.size(); ++downspin_idx)
	    {
	      uint64 send_idx = n_upspins_i_send_back_offsets_[destination_mpi_rank] +
		n_states_already_prepared[destination_mpi_rank]++; 
	      uint64 downspin_offset = my_downspins_offset_[my_downspins_[downspin_idx]];
	      uint64 idx = upspin_offset + downspin_offset;
	      send_buffer_[send_idx] = recv_buffer_[idx];
	    }
	  ++upspin_offset;
	}
      t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  send back (prepare):   %3.4f\n", t2-t1); 




      // // DEBUG START
      // MPI_Barrier(MPI_COMM_WORLD);
      // if(mpi_rank_ == 0) printf("BACK SEND BUFFER --------------------------------\n");
      // MPI_Barrier(MPI_COMM_WORLD);
      // for(int nt = 0; nt < mpi_size_; ++nt) 
      //   {
      //     if(mpi_rank_ == nt) 
      //       {
      // 	printf("[%d]:\n", mpi_rank_);

      // 	for (int sendto = 0; sendto < mpi_size_; ++sendto)
      // 	  {
      // 	    int from = n_upspins_i_send_back_offsets_[sendto];
      // 	    int to = (sendto == mpi_size_-1) ? sum_n_upspins_i_send_back_ : n_upspins_i_send_back_offsets_[sendto+1];
      // 	    for (int idx = from; idx < to; ++idx)
      // 	      printf("%f -> [%d] %s; %s \n", 
      // 		     send_buffer_[idx], sendto,
      // 		     PrintSpinhalf(n_sites_, upspins_i_send_back_[idx]).c_str(),
      // 		     PrintSpinhalf(n_sites_, downspins_i_send_back_[idx]).c_str());
      // 	  }
      // 	printf("\n");
      //       }
      //     MPI_Barrier(MPI_COMM_WORLD);
      //   }
      // MPI_Barrier(MPI_COMM_WORLD);
      // // DEBUG END


      t1 = MPI_Wtime();
      lila::MPI_Alltoallv<coeff_t>
	(send_buffer_.data(), n_upspins_i_send_back_.data(), 
	 n_upspins_i_send_back_offsets_.data(), 
	 recv_buffer_.data(), n_upspins_i_recv_back_.data(), 
	 n_upspins_i_recv_back_offsets_.data(), 
	 MPI_COMM_WORLD);
      t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  send back (Alltoall):   %3.4f\n", t2-t1); 



      // // DEBUG START
      // MPI_Barrier(MPI_COMM_WORLD);
      // if(mpi_rank_ == 0) printf("BACK RECV BUFFER --------------------------------\n");
      // MPI_Barrier(MPI_COMM_WORLD);
      // for(int nt = 0; nt < mpi_size_; ++nt) 
      //   {
      //     if(mpi_rank_ == nt) 
      //       {
      // 	printf("[%d]:\n", mpi_rank_);
      // 	for (int recvfrom = 0; recvfrom < mpi_size_; ++recvfrom)
      // 	  {
      // 	    int from = n_upspins_i_recv_back_offsets_[recvfrom];
      // 	    int to = (recvfrom == mpi_size_-1) ? sum_n_upspins_i_recv_back_ : n_upspins_i_recv_back_offsets_[recvfrom+1];
      // 	    for (int idx = from; idx < to; ++idx)
      // 	      printf("%f <- [%d]  %s; %s\n", recv_buffer_[idx], recvfrom,
      // 		     PrintSpinhalf(n_sites_, upspins_i_recv_back_[idx]).c_str(),
      // 		     PrintSpinhalf(n_sites_, downspins_i_recv_back_[idx]).c_str());
		    
      // 	  }
      // 	printf("\n");
      //       }
      //     MPI_Barrier(MPI_COMM_WORLD);
      //   }
      // MPI_Barrier(MPI_COMM_WORLD);
      // // DEBUG END



      // Sort the received data into the send buffer
      t1 = MPI_Wtime();
      std::fill(send_buffer_.begin(), send_buffer_.end(), 0);
      for (uint64 idx = 0; idx < sum_n_upspins_i_recv_back_; ++idx)
	{
	  state_t upspins = upspins_i_recv_back_[idx];
	  state_t downspins = downspins_i_recv_back_[idx];	    
	  uint64 sorted_idx = my_upspins_offset_[upspins] + 
	    indexing_downspins_.index(downspins);
	  send_buffer_[sorted_idx] = recv_buffer_[idx];
	}

      t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  sort back:   %3.4f\n", t2-t1); 



      // // DEBUG START
      // MPI_Barrier(MPI_COMM_WORLD);
      // if(mpi_rank_ == 0) printf("BACK SORTED SEND BUFFER --------------------------------\n");
      // MPI_Barrier(MPI_COMM_WORLD);
      // for(int nt = 0; nt < mpi_size_; ++nt) 
      //   {
      //     if(mpi_rank_ == nt) 
      //       {
      // 	printf("[%d]:\n", mpi_rank_);
      // 	uint64 idx = 0;
      // 	for (const state_t& upspins : my_upspins_) 
      // 	  {
      // 	    for (state_t downspins : hs_downspins_)
      // 	      {
      // 		printf("%s; %s  %f\n", 
      // 		       PrintSpinhalf(n_sites_, upspins).c_str(),
      // 		       PrintSpinhalf(n_sites_, downspins).c_str(),
      // 		       send_buffer_[idx]);
      // 		++idx;
      // 	      }
      // 	  }
      // 	printf("\n");
      //       }
      //     MPI_Barrier(MPI_COMM_WORLD);
      //   }
      // MPI_Barrier(MPI_COMM_WORLD);
      // // DEBUG END


      // // DEBUG START
      // // DEBUG FORWARD_BACKWARD COMMUNICATION
      // for (int k=0; k<in_vec.vector_local().size(); ++k)
      //   assert(send_buffer_[k] == in_vec.vector_local()(k));
      // printf("[%d] SUCCESS !!!!!!!!!!!!!!!!!\n", mpi_rank_);
      // MPI_Finalize();
      // exit(1);
      // // DEBUG END


      // Fill in from send buffer
      t1 = MPI_Wtime();
      for (int k=0; k<in_vec.vector_local().size(); ++k)
	out_vec.vector_local()(k) += send_buffer_[k];
      t2 = MPI_Wtime();
      if ((mpi_rank_ == 0) && verbose) printf("  fill outvec:   %3.4f\n", t2-t1); 

    }

    template <class coeff_t, class state_t>
    hubbard_qn HubbardModelMPI<coeff_t, state_t>::
    apply_fermion
    (const lila::VectorMPI<coeff_t>& in_vec, lila::VectorMPI<coeff_t>& out_vec, 
     std::string type, int site)
    {
      using utils::popcnt;
      using utils::gbit;
      using utils::gbits;
      assert(site < n_sites_);

      hubbard_qn qn_after = qn_;
      if (type == "cdagup") 
	{
	  ++qn_after.n_upspins;
	  if (mpi_rank_== 0) 
	    printf("cdagup not implemented! Use cdagdn instead.\n");
	  MPI_Finalize();
	  exit(EXIT_FAILURE);
	}
      else if (type == "cup")
	{
	  --qn_after.n_upspins;
	  if (mpi_rank_== 0) 
	    printf("cup not implemented! Use cdn instead.\n");
	  MPI_Finalize();
	  exit(EXIT_FAILURE);
	}
      else if (type == "cdagdn") ++qn_after.n_downspins;
      else if (type == "cdn") --qn_after.n_downspins;
      else 
	{
	  if (mpi_rank_== 0) 
	    printf("Error in apply_fermion: Invalid fermion type!\n");
	  MPI_Finalize();
	  exit(EXIT_FAILURE);
	}

      Spinhalf<state_t> hs_downspins_after(n_sites_, qn_after.n_downspins);
      IndexTable<Spinhalf<state_t>, uint64> 
	indexing_downspins_after(hs_downspins_after);

      // Compute new offsets for altered hilbertspace
      std::unordered_map<state_t, uint64> my_upspins_offset_after;
      for (auto old_offset : my_upspins_offset_)
	{
	  state_t state = old_offset.first;
	  uint64 off = old_offset.second;
	  uint64 new_off = off / indexing_downspins_.size() 
	    * indexing_downspins_after.size();
	  my_upspins_offset_after[state] = new_off;
	}

      // Allocate out_vec
      uint64 local_dim_after = 
	my_upspins_.size() * indexing_downspins_after.size();
      try 
	{
	  out_vec.resize(local_dim_after);
	}
      catch(...)
	{
	  std::cerr << "[ " << mpi_rank_
		    << " ] Error: Could not allocate fermion vector!" 
		    << std::endl << std::flush;
	  MPI_Abort(MPI_COMM_WORLD, 4);
	}   


      // Loop over all configurations
      uint64 upspin_idx = 0;
      for (const state_t& upspins : my_upspins_) 
	{
	  uint64 upspin_offset = my_upspins_offset_[upspins];
	  uint64 upspin_offset_after = my_upspins_offset_after[upspins];
	  uint64 downspin_offset = 0;


	  // Apply the electron part
	  if (type == "cdagdn")
	    {
	      const state_t sitemask = ((state_t)1 << site);

	      for (state_t downspins : hs_downspins_)
		{

		  // raise local site val if 0
		  if (gbit(downspins, site) == 0) 
		    {
		      uint64 idx = upspin_offset + downspin_offset;
		      state_t new_downspins = downspins;
		      new_downspins |= sitemask;
		      uint64 new_idx = upspin_offset_after
			+ indexing_downspins_after.index(new_downspins);
		      double fermi = 
			popcnt(gbits(downspins, site, 0)) % 2==0 ? 1. : -1.;
		      out_vec.vector_local()(new_idx) += fermi * in_vec.vector_local()(idx);
		    }

		  ++downspin_offset;
		}
	    }  // cdagdn


	  // Apply the hole part
	  else if (type == "cdn")
	    {
	      const state_t antisitemask = ~((state_t)1 << site);

	      for (state_t downspins : hs_downspins_)
		{

		  // lower local site val if 1
		  if (gbit(downspins, site) == 1) 
		    {
		      uint64 idx = upspin_offset + downspin_offset;
		      state_t new_downspins = downspins;
		      new_downspins &= antisitemask;
		      uint64 new_idx = upspin_offset_after
			+ indexing_downspins_after.index(new_downspins);
		      double fermi = 
			popcnt(gbits(downspins, site, 0)) % 2==0 ? 1. : -1.;
		      out_vec.vector_local()(new_idx) += fermi * in_vec.vector_local()(idx);
		    }

		  ++downspin_offset;
		}
	    }  // cdn


	  ++upspin_idx;
	}  // loop over upspin configurations

      return qn_after;
    }

      
    template <class coeff_t, class state_t>
    void HubbardModelMPI<coeff_t, state_t>::initialize()
    {

      hs_upspins_ = Spinhalf<state_t>(n_sites_, qn_.n_upspins);
      hs_downspins_ = Spinhalf<state_t>(n_sites_, qn_.n_downspins);
      
      indexing_upspins_= IndexTable<Spinhalf<state_t>, uint64>(hs_upspins_);
      indexing_downspins_ = IndexTable<Spinhalf<state_t>, uint64>(hs_downspins_);

      // clear previous members
      my_upspins_.clear();
      my_upspins_offset_.clear();
      my_downspins_.clear();
      my_downspins_offset_.clear();
      n_downspins_i_send_forward_.clear();
      n_downspins_i_recv_forward_.clear();
      n_downspins_i_send_forward_offsets_.clear();
      n_downspins_i_recv_forward_offsets_.clear();

      downspins_i_recv_forward_.clear();
      upspins_i_recv_forward_.clear();
      downspins_i_send_forward_offsets_.clear();
      downspins_i_recv_forward_offsets_.clear();

      sum_n_downspins_i_send_forward_ = 0;
      sum_n_downspins_i_recv_forward_ = 0; 

      n_upspins_i_send_back_.clear();
      n_upspins_i_recv_back_.clear();
      n_upspins_i_send_back_offsets_.clear();
      n_upspins_i_recv_back_offsets_.clear();
      downspins_i_recv_back_.clear();
      upspins_i_recv_back_.clear();
      upspins_i_send_back_offsets_.clear();
      upspins_i_recv_back_offsets_.clear();

      sum_n_upspins_i_send_back_ = 0;
      sum_n_upspins_i_recv_back_ = 0; 

      buffer_size_ = 0;
      send_buffer_.clear();
      recv_buffer_.clear();
 
      // Collect all upspin configurations belonging to this mpi_rank
      uint64 offset = 0;
      for (auto state : hs_upspins_)
	{
	  if (mpi_rank_of_spins(state) == mpi_rank_) 
	    {
	      my_upspins_.push_back(state);
	      my_upspins_offset_[state] = offset;	
	      offset += hs_downspins_.size();
	    }
	}

      uint64 dim_test = my_upspins_.size();
      MPI_Allreduce(MPI_IN_PLACE, &dim_test, 1, MPI_UNSIGNED_LONG, MPI_SUM, 
		    MPI_COMM_WORLD);
      assert(dim_test == hs_upspins_.size());
	
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
	      offset += hs_upspins_.size();
	    }
	}
      local_dim_downspins_ = offset;
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

      for (state_t downspins : hs_downspins_)
	{
	  int destination_mpi_rank = mpi_rank_of_spins(downspins);
	  n_downspins_i_send_forward_[destination_mpi_rank] += my_upspins_.size();
	}
      // Communicate number of send/receive states
      MPI_Alltoall(n_downspins_i_send_forward_.data(), 1, MPI_INT,
		   n_downspins_i_recv_forward_.data(), 1, MPI_INT, MPI_COMM_WORLD);


      // Compute offsets and total number of send/receive states
      for(int m = 0; m < mpi_size_; ++m)
	for(int n = 0; n < m; ++n)
	  {
	    n_downspins_i_send_forward_offsets_[m] += n_downspins_i_send_forward_[n];
	    n_downspins_i_recv_forward_offsets_[m] += n_downspins_i_recv_forward_[n];
	  }

      sum_n_downspins_i_send_forward_ = std::accumulate(n_downspins_i_send_forward_.begin(),
							n_downspins_i_send_forward_.end(), 0);
      sum_n_downspins_i_recv_forward_ = std::accumulate(n_downspins_i_recv_forward_.begin(),
							n_downspins_i_recv_forward_.end(), 0);

      // Communicate which downstates are sent/received by whom
      std::vector<state_t> downspins_i_send_forward_(sum_n_downspins_i_send_forward_, 0);
      std::vector<state_t> upspins_i_send_forward_(sum_n_downspins_i_send_forward_, 0);

      // downspins_i_send_forward_ = std::vector<state_t>(sum_n_downspins_i_send_forward_, 0);
      // upspins_i_send_forward_= std::vector<state_t>(sum_n_downspins_i_send_forward_, 0);

      std::vector<uint64> n_downspins_already_prepared(mpi_size_, 0);
      for (state_t downspins : hs_downspins_)
	{	
	  int destination_mpi_rank = mpi_rank_of_spins(downspins);
	  for (const state_t& upspins : my_upspins_)  // loop over upspins of process
	    {
	      uint64 idx = n_downspins_i_send_forward_offsets_[destination_mpi_rank] +
		n_downspins_already_prepared[destination_mpi_rank]++; 
	      downspins_i_send_forward_[idx] = downspins;
	      upspins_i_send_forward_[idx] = upspins;
	    }
	}
      downspins_i_recv_forward_.resize(sum_n_downspins_i_recv_forward_);
      upspins_i_recv_forward_.resize(sum_n_downspins_i_recv_forward_);

      MPI_Alltoallv(downspins_i_send_forward_.data(), 
		    n_downspins_i_send_forward_.data(), n_downspins_i_send_forward_offsets_.data(), 
		    MPI_UNSIGNED,
		    downspins_i_recv_forward_.data(), 		      
		    n_downspins_i_recv_forward_.data(), n_downspins_i_recv_forward_offsets_.data(), 
		    MPI_UNSIGNED, MPI_COMM_WORLD);
      MPI_Alltoallv(upspins_i_send_forward_.data(), 
		    n_downspins_i_send_forward_.data(), n_downspins_i_send_forward_offsets_.data(), 
		    MPI_UNSIGNED,
		    upspins_i_recv_forward_.data(), 		      
		    n_downspins_i_recv_forward_.data(), n_downspins_i_recv_forward_offsets_.data(), 
		    MPI_UNSIGNED, MPI_COMM_WORLD);


      // downspins_i_send_forward_ and upspins_i_send_forward_ not neede anymore
      downspins_i_send_forward_.clear();
      upspins_i_send_forward_.clear();
      downspins_i_send_forward_.shrink_to_fit();
      upspins_i_send_forward_.shrink_to_fit();

      // // Print forward communication patterns
      // MPI_Barrier(MPI_COMM_WORLD);
      // for(int nt = 0; nt < mpi_size_; ++nt) 
      //   {
      //     if(mpi_rank_ == nt) 
      //       {
      // 	printf("[%d]: forward\n", mpi_rank_);
      // 	printf("my_upspins.size(): %d\n", my_upspins_.size());				      
      // 	printf("sum_n_downspins_i_send_forward_: %d\n", sum_n_downspins_i_send_forward_);				      
      // 	printf("sum_n_downspins_i_recv_forward_: %d\n", sum_n_downspins_i_recv_forward_);
	
      // 	for(int k = 0; k < mpi_size_; ++k)
      // 	  { 
      // 	    printf("-> %d (down): %d\n", k, n_downspins_i_send_forward_[k]);
      // 	    int from = n_downspins_i_send_forward_offsets_[k];
      // 	    int to = (k == mpi_size_-1) ? sum_n_downspins_i_send_forward_ : n_downspins_i_send_forward_offsets_[k+1];
      // 	    for (int l = from; l < to; ++l)
      // 	      printf("  %s \n", PrintSpinhalf(n_sites_, downspins_i_send_forward_[l]).c_str());
      // 	  }

      // 	for(int k = 0; k < mpi_size_; ++k) 
      // 	  {
      // 	    printf("<- %d (down): %d\n", k, n_downspins_i_recv_forward_[k]);
      // 	    int from = n_downspins_i_recv_forward_offsets_[k];
      // 	    int to = (k == mpi_size_-1) ? sum_n_downspins_i_recv_forward_ : n_downspins_i_recv_forward_offsets_[k+1];
      // 	    for (int l = from; l < to; ++l)
      // 	      printf("  %s \n", PrintSpinhalf(n_sites_, downspins_i_recv_forward_[l]).c_str());		  
      // 	  }
      // 	printf("\n");
      //       }
      //     MPI_Barrier(MPI_COMM_WORLD);
      //   }
      // MPI_Barrier(MPI_COMM_WORLD);




      ///////////////////////////////////////////////
      // Determine Communication patterns
      // for back sending downspins
      n_upspins_i_send_back_.resize(mpi_size_, 0);
      n_upspins_i_recv_back_.resize(mpi_size_, 0);
      n_upspins_i_send_back_offsets_.resize(mpi_size_, 0);
      n_upspins_i_recv_back_offsets_.resize(mpi_size_, 0);

      for (state_t upspins : hs_upspins_)
	{
	  int destination_mpi_rank = mpi_rank_of_spins(upspins);
	  n_upspins_i_send_back_[destination_mpi_rank] += my_downspins_.size();
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

      // downspins_i_send_back_ = std::vector<state_t>(sum_n_upspins_i_send_back_, 0);
      // upspins_i_send_back_ = std::vector<state_t>(sum_n_upspins_i_send_back_, 0);

      std::vector<uint64> n_upspins_already_prepared(mpi_size_, 0);
      for (state_t upspins : hs_upspins_)
	{	
	  int destination_mpi_rank = mpi_rank_of_spins(upspins);
	  for (const state_t& downspins : my_downspins_)  // loop over upspins of process
	    {
	      uint64 idx = n_upspins_i_send_back_offsets_[destination_mpi_rank] +
		n_upspins_already_prepared[destination_mpi_rank]++; 
	      downspins_i_send_back_[idx] = downspins;
	      upspins_i_send_back_[idx] = upspins;
	    }
	}
      downspins_i_recv_back_.resize(sum_n_upspins_i_recv_back_);
      upspins_i_recv_back_.resize(sum_n_upspins_i_recv_back_);

      MPI_Alltoallv(downspins_i_send_back_.data(), 
		    n_upspins_i_send_back_.data(), n_upspins_i_send_back_offsets_.data(), 
		    MPI_UNSIGNED,
		    downspins_i_recv_back_.data(), 		      
		    n_upspins_i_recv_back_.data(), n_upspins_i_recv_back_offsets_.data(), 
		    MPI_UNSIGNED, MPI_COMM_WORLD);
      MPI_Alltoallv(upspins_i_send_back_.data(), 
		    n_upspins_i_send_back_.data(), n_upspins_i_send_back_offsets_.data(), 
		    MPI_UNSIGNED,
		    upspins_i_recv_back_.data(), 		      
		    n_upspins_i_recv_back_.data(), n_upspins_i_recv_back_offsets_.data(), 
		    MPI_UNSIGNED, MPI_COMM_WORLD);

      // downspins_i_send_back_ and upspins_i_send_back_ not neede anymore
      downspins_i_send_back_.clear();
      upspins_i_send_back_.clear();
      downspins_i_send_back_.shrink_to_fit();
      upspins_i_send_back_.shrink_to_fit();


      // // Print back communication patterns
      // MPI_Barrier(MPI_COMM_WORLD);
      // for(int nt = 0; nt < mpi_size_; ++nt) 
      //   {
      //     if(mpi_rank_ == nt) 
      //       {
      // 	printf("[%d]: back\n", mpi_rank_);
      // 	printf("my_upspins.size(): %d\n", my_upspins_.size());				      
      // 	printf("sum_n_upspins_i_send_back_: %d\n", sum_n_upspins_i_send_back_);				      
      // 	printf("sum_n_upspins_i_recv_back_: %d\n", sum_n_upspins_i_recv_back_);
	
      // 	for(int k = 0; k < mpi_size_; ++k)
      // 	  { 
      // 	    printf("-> %d (up): %d\n", k, n_upspins_i_send_back_[k]);
      // 	    int from = n_upspins_i_send_back_offsets_[k];
      // 	    int to = (k == mpi_size_-1) ? sum_n_upspins_i_send_back_ : n_upspins_i_send_back_offsets_[k+1];
      // 	    for (int l = from; l < to; ++l)
      // 	      printf("  %s; %s\n", PrintSpinhalf(n_sites_, upspins_i_send_back_[l]).c_str(),
      // 		     PrintSpinhalf(n_sites_, downspins_i_send_back_[l]).c_str());
      // 	  }

      // 	for(int k = 0; k < mpi_size_; ++k) 
      // 	  {
      // 	    printf("<- %d (up): %d\n", k, n_upspins_i_recv_back_[k]);
      // 	    int from = n_upspins_i_recv_back_offsets_[k];
      // 	    int to = (k == mpi_size_-1) ? sum_n_upspins_i_recv_back_ : n_upspins_i_recv_back_offsets_[k+1];
      // 	    for (int l = from; l < to; ++l)
      // 	      printf("  %s \n", PrintSpinhalf(n_sites_, upspins_i_recv_back_[l]).c_str());		  
      // 	  }
      // 	printf("\n");
      //       }
      //     MPI_Barrier(MPI_COMM_WORLD);
      //   }
      // MPI_Barrier(MPI_COMM_WORLD);


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


    template class HubbardModelMPI<double, uint32>;
    template class HubbardModelMPI<complex, uint32>;
    template class HubbardModelMPI<double, uint64>;
    template class HubbardModelMPI<complex, uint64>;
       
  }  // namespace models
}  // namespace hydra
