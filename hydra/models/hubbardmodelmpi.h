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

// #include <functional>
// #include <utility>
#include <vector>
#include <unordered_map>
#include <lila/allmpi.h>

#include <hydra/hilbertspaces/hubbard.h>
#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/indexing/indextable.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra { namespace models {
    
    using hilbertspaces::Spinhalf;
    using hilbertspaces::PrintSpinhalf;
    using hilbertspaces::hubbard_qn;
    using indexing::IndexTable;
    using operators::BondList;
    using operators::Couplings;

    template <class coeff_t, class state_t=uint32>
    class HubbardModelMPI
    {
    public:
      /*! 
	Defines a Hubbard model given number of sites and pairs of 
	neighboring sites with hoppings and initializes communication
      */
      HubbardModelMPI(BondList bondlist, Couplings couplings, 
		      hilbertspaces::hubbard_qn qn);
      
      void apply_hamiltonian(const lila::VectorMPI<coeff_t>& in_vec,
			     lila::VectorMPI<coeff_t>& out_vec,
			     bool verbose = false);

      hilbertspaces::hubbard_qn apply_fermion
      (const lila::VectorMPI<coeff_t>& in_vec, lila::VectorMPI<coeff_t>& out_vec,
       std::string type, int site);

      hilbertspaces::hubbard_qn qn() const { return qn_; }
      void set_qn(hilbertspaces::hubbard_qn qn);
      int n_sites() const { return n_sites_; }
      uint64 local_dim() const { return local_dim_; }
      uint64 dim() const { return dim_; }

    private:     

      const int n_sites_;
      hubbard_qn qn_;

      std::vector<std::pair<int, int>> hoppings_;
      std::vector<coeff_t> hopping_amplitudes_;
      std::vector<std::pair<int, int>> currents_;
      std::vector<coeff_t> current_amplitudes_;
      std::vector<std::pair<int, int>> interactions_;
      std::vector<double> interaction_strengths_;
      std::vector<int> onsites_;
      std::vector<double> onsite_potentials_;
      double U_;

      // members below are all set by initialize
      void initialize();
      Spinhalf<state_t> hs_upspins_;
      Spinhalf<state_t> hs_downspins_;
      IndexTable<Spinhalf<state_t>, uint64> indexing_upspins_;
      IndexTable<Spinhalf<state_t>, uint64> indexing_downspins_;

      int mpi_rank_;
      int mpi_size_;

      uint64 dim_;
      uint64 local_dim_;
      uint64 max_local_dim_;
      uint64 min_local_dim_;
      uint64 local_dim_downspins_;

      std::vector<state_t> my_upspins_;
      std::unordered_map<state_t, uint64> my_upspins_offset_;
      std::vector<state_t> my_downspins_;
      std::unordered_map<state_t, uint64> my_downspins_offset_;
      
      // Determine process of a spin configuration
      inline int mpi_rank_of_spins(const state_t& spins) const
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
      std::vector<coeff_t> send_buffer_;
      std::vector<coeff_t> recv_buffer_;

      // std::vector<state_t> downspins_i_send_forward_;
      // std::vector<state_t> upspins_i_send_forward_;
      // std::vector<state_t> downspins_i_send_back_;
      // std::vector<state_t> upspins_i_send_back_;
     

    };
    
    




  }
}


#endif
