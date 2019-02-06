// Copyright 2019 Alexander Wietek - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef HYDRA_MODELS_HUBBARDMODELMPI_
#define HYDRA_MODELS_HUBBARDMODELMPI_

#include <functional>
#include <utility>
#include <vector>
#include <unordered_map>

#include <lila/all.h>
#include <hydra/hilbertspaces/hubbard.h>

namespace hydra { namespace models {
    
    using hilbertspaces::Spinhalf;
    using hilbertspaces::PrintSpinhalf;
    using hilbertspaces::hubbard_qn;
    using indexing::IndexTable;

    template <class state_t>
    class HubbardModelMPI
    {
    public:
      HubbardModelMPI(const int& n_sites, 
		      const std::vector<std::pair<int, int>> neighbors,
		      const hubbard_qn& qn)
	: n_sites_(n_sites),
	  neighbors_(neighbors),
	  qn_(qn),
	  hs_upspins_(n_sites_, qn_.n_upspins),
	  hs_downspins_(n_sites_, qn_.n_downspins),
	  indexing_upspins_(hs_upspins_),
	  indexing_downspins_(hs_downspins_)
      {
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
	// if (mpi_rank_ == 0) printf("Sleeping constructor model I ...\n");
	// sleep(20);

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

	// if (mpi_rank_ == 0) printf("Sleeping constructor model II ...\n");
	// sleep(20);

	uint64 dim_test = my_upspins_.size();
	MPI_Allreduce(MPI_IN_PLACE, &dim_test, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	assert(dim_test == hs_upspins_.size());
	
	local_dim_ = offset;
	MPI_Allreduce(&local_dim_, &dim_, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	

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


	// if (mpi_rank_ == 0) printf("Sleeping constructor model III ...\n");
	// sleep(20);

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

	// if (mpi_rank_ == 0) printf("Sleeping constructor model IV ...\n");
	// sleep(20);

	// Communicate which downstates are sent/received by whom
	std::vector<state_t> downspins_i_send_forward_(sum_n_downspins_i_send_forward_, 0);
	std::vector<state_t> upspins_i_send_forward_(sum_n_downspins_i_send_forward_, 0);

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

	// if (mpi_rank_ == 0) printf("Sleeping constructor model Va ...\n");
	// sleep(20);

	// downspins_i_send_forward_ and upspins_i_send_forward_ not neede anymore
	downspins_i_send_forward_.clear();
	upspins_i_send_forward_.clear();
	downspins_i_send_forward_.shrink_to_fit();
	upspins_i_send_forward_.shrink_to_fit();

	// if (mpi_rank_ == 0) printf("Sleeping constructor model V ...\n");
	// sleep(20);

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

	// if (mpi_rank_ == 0) printf("Sleeping constructor model VI ...\n");
	// sleep(20);

	// Communicate which upstates are sent/received by whom

	std::vector<state_t> downspins_i_send_back_(sum_n_upspins_i_send_back_, 0);
	std::vector<state_t> upspins_i_send_back_(sum_n_upspins_i_send_back_, 0);

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

	// if (mpi_rank_ == 0) printf("Sleeping constructor model VIIa ...\n");
	// sleep(20);

	// downspins_i_send_back_ and upspins_i_send_back_ not neede anymore
	downspins_i_send_back_.clear();
	upspins_i_send_back_.clear();
	downspins_i_send_back_.shrink_to_fit();
	upspins_i_send_back_.shrink_to_fit();

	// if (mpi_rank_ == 0) printf("Sleeping constructor model VII ...\n");
	// sleep(20);

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


	// if (mpi_rank_ == 0) printf("Sleeping constructor model VIII ...\n");
	// sleep(20);
 

      }
      
      void apply_hamiltonian(const double& t, const double& U,
			     const lila::VectorMPI<double>& in_vec,
			     lila::VectorMPI<double>& out_vec)
      {
	using utils::popcnt;
	using utils::gbits;

	auto hs_downspins = Spinhalf<state_t>(n_sites_, qn_.n_downspins);

	Zeros(out_vec.vector_local());


	double t1 = MPI_Wtime();
	// Apply Diagonal terms
	uint64 upspin_idx = 0;
	for (const state_t& upspins : my_upspins_)  // loop over upspins of process
	  {
	    uint64 upspin_offset = my_upspins_offset_[upspins];
	    uint64 downspin_offset=0;
	    for (state_t downspins : hs_downspins)  // loop over all downspins
	      {
		uint64 idx = upspin_offset + downspin_offset;
		double coeff = U*(double)popcnt(upspins & downspins);
		out_vec.vector_local()(idx) += coeff*in_vec.vector_local()(idx);
		++downspin_offset;
	      }
	    ++upspin_idx;
	  }
	double t2 = MPI_Wtime();
	if (mpi_rank_== 0) printf("  diag: %3.4f\n", t2-t1); 


      
	// Apply hoppings on downspins
	t1 = MPI_Wtime();
	for (auto pair : neighbors_)
	  {
	    int s1 = std::min(pair.first, pair.second); 
	    int s2 = std::max(pair.first, pair.second);
	    uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);

	    // Loop over all configurations
	    upspin_idx = 0;
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
	  }  // hopping on downspin
	t2 = MPI_Wtime();
	if (mpi_rank_== 0) printf("  down: %3.4f\n", t2-t1); 




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
	    for (upspin_idx = 0; upspin_idx < my_upspins_.size(); ++upspin_idx)
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
	if (mpi_rank_== 0) printf("  send forward (prepare): %3.4f\n", t2-t1); 
	
	t1 = MPI_Wtime();
	MPI_Alltoallv(send_buffer_.data(), 
		      n_downspins_i_send_forward_.data(), 
		      n_downspins_i_send_forward_offsets_.data(), 
		      MPI_DOUBLE,
		      recv_buffer_.data(),
		      n_downspins_i_recv_forward_.data(), 
		      n_downspins_i_recv_forward_offsets_.data(), 
		      MPI_DOUBLE, MPI_COMM_WORLD);
	t2 = MPI_Wtime();
	if (mpi_rank_== 0) printf("  send forward (Alltoallv): %3.4f\n", t2-t1); 
	
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
	if (mpi_rank_== 0) printf("  sort forward: %3.4f\n", t2-t1); 

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

	
	// Apply upspin hoppings on the sorted send buffer
	// write the results to the receive buffer
	
	// // JUST COPY INVEC FOR DEBUGGING
	// std::copy(send_buffer_.begin(), send_buffer_.end(), recv_buffer_.begin());
	

	// Apply hoppings on upspins
	t1 = MPI_Wtime();
	std::fill(recv_buffer_.begin(), recv_buffer_.end(), 0);
	for (auto pair : neighbors_)
	  {
	    int s1 = std::min(pair.first, pair.second); 
	    int s2 = std::max(pair.first, pair.second);
	    uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);

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
	  }  // hopping on downspin
	t2 = MPI_Wtime();
	if (mpi_rank_== 0) printf("  up:   %3.4f\n", t2-t1); 
	


	
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
	if (mpi_rank_== 0) printf("  send back (prepare):   %3.4f\n", t2-t1); 

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

	t1 = MPI_Wtime();
	MPI_Alltoallv(send_buffer_.data(), 
		      n_downspins_i_send_forward_.data(), 
		      n_downspins_i_send_forward_offsets_.data(), 
		      MPI_DOUBLE,
		      recv_buffer_.data(),
		      n_downspins_i_recv_forward_.data(), 
		      n_downspins_i_recv_forward_offsets_.data(), 
		      MPI_DOUBLE, MPI_COMM_WORLD);
	t2 = MPI_Wtime();
	if (mpi_rank_== 0) printf("  send back (Alltoall):   %3.4f\n", t2-t1); 

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


	// Sort the received data into the send buffer
	t1 = MPI_Wtime();
	std::fill(send_buffer_.begin(), send_buffer_.end(), 0);
	for (uint64 idx = 0; idx < sum_n_downspins_i_recv_forward_; ++idx)
	  {
	    state_t upspins = upspins_i_recv_back_[idx];
	    state_t downspins = downspins_i_recv_back_[idx];	    
	    
	    uint64 sorted_idx = my_upspins_offset_[upspins] + 
	      indexing_downspins_.index(downspins);
	    send_buffer_[sorted_idx] = recv_buffer_[idx];
	  }
	t2 = MPI_Wtime();
	if (mpi_rank_== 0) printf("  sort back:   %3.4f\n", t2-t1); 

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


	// // DEBUG FORWARD_BACKWARD COMMUNICATION
	// for (int k=0; k<in_vec.vector_local().size(); ++k)
	//   assert(send_buffer_[k] == in_vec.vector_local()(k));
	// printf("[%d] SUCCESS !!!!!!!!!!!!!!!!!\n", mpi_rank_);

	// MPI_Finalize();
	// exit(1);

	// Fill in from send buffer
	t1 = MPI_Wtime();
	for (int k=0; k<in_vec.vector_local().size(); ++k)
	  out_vec.vector_local()(k) += send_buffer_[k];
	t2 = MPI_Wtime();
	if (mpi_rank_== 0) printf("  fill outvec:   %3.4f\n", t2-t1); 

      }
      
      uint64 local_dim() const { return local_dim_; }
      uint64 dim() const { return dim_; }


    private:     
      const int n_sites_;
      std::vector<std::pair<int, int>> neighbors_;
      hubbard_qn qn_;
      Spinhalf<state_t> hs_upspins_;
      Spinhalf<state_t> hs_downspins_;
      IndexTable<Spinhalf<state_t>, uint64> indexing_upspins_;
      IndexTable<Spinhalf<state_t>, uint64> indexing_downspins_;

      int mpi_rank_;
      int mpi_size_;

      uint64 dim_;
      uint64 local_dim_;
      uint64 local_dim_downspins_;

      std::vector<state_t> my_upspins_;
      std::unordered_map<state_t, uint64> my_upspins_offset_;
      std::vector<state_t> my_downspins_;
      std::unordered_map<state_t, uint64> my_downspins_offset_;
      
      // Determine process of a spin configuration
      int mpi_rank_of_spins(const state_t& spins) const
      { 
	// return (int)(std::hash<state_t>{}(spins) % mpi_size_); 
	state_t x = spins;
	x = ((x >> 16) ^ x) * 0x45d9f3b;
	x = ((x >> 16) ^ x) * 0x45d9f3b;
	x = (x >> 16) ^ x;
	return x % mpi_size_;
      }

      
      // Communication patterns
      std::vector<int> n_downspins_i_send_forward_;
      std::vector<int> n_downspins_i_recv_forward_;
      std::vector<int> n_downspins_i_send_forward_offsets_;
      std::vector<int> n_downspins_i_recv_forward_offsets_;
      std::vector<state_t> downspins_i_recv_forward_;
      std::vector<state_t> upspins_i_recv_forward_;
      std::vector<int> downspins_i_send_forward_offsets_;
      std::vector<int> downspins_i_recv_forward_offsets_;

      uint64 sum_n_downspins_i_send_forward_;
      uint64 sum_n_downspins_i_recv_forward_; 

      std::vector<int> n_upspins_i_send_back_;
      std::vector<int> n_upspins_i_recv_back_;
      std::vector<int> n_upspins_i_send_back_offsets_;
      std::vector<int> n_upspins_i_recv_back_offsets_;
      std::vector<state_t> downspins_i_recv_back_;
      std::vector<state_t> upspins_i_recv_back_;
      std::vector<int> upspins_i_send_back_offsets_;
      std::vector<int> upspins_i_recv_back_offsets_;

      uint64 sum_n_upspins_i_send_back_;
      uint64 sum_n_upspins_i_recv_back_; 

      uint64 buffer_size_;
      std::vector<double> send_buffer_;
      std::vector<double> recv_buffer_;

    };
    
  }
}


#endif
