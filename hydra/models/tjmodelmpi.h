// Copyright 2020 Alexander Wietek - All Rights Reserved.
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

#ifndef HYDRA_MODELS_TJMODELMPI_
#define HYDRA_MODELS_TJMODELMPI_

#include <vector>
#include <unordered_map>

#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/hilbertspaces/hubbard.h>
#include <hydra/indexing/indextable.h>
#include <hydra/indexing/lintable.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <lila/allmpi.h>

using namespace hydra::hilbertspaces;
using namespace hydra::indexing;
using namespace hydra::operators;

namespace hydra { namespace models {


    template <class coeff_t, class state_t=uint32>
    class TJModelMPI
    {
    public:

      TJModelMPI(BondList bondlist, Couplings couplings, 
		 hilbertspaces::hubbard_qn qn);
      
      void apply_hamiltonian(lila::VectorMPI<coeff_t> const& in_vec,
			     lila::VectorMPI<coeff_t>& out_vec,
			     bool verbose = false);

      hubbard_qn qn() const { return qn_; }
      int n_sites() const { return n_sites_; }
      uint64 local_dim() const { return local_dim_; }
      uint64 dim() const { return dim_; }

    private:     
      const int n_sites_;
      hubbard_qn qn_;

      int mpi_rank_;
      int mpi_size_;

      // Dimensions
      uint64 dim_;
      uint64 local_dim_;
      uint64 max_local_dim_;
      uint64 min_local_dim_;
      uint64 local_dim_downspins_;

      std::vector<std::pair<int, int>> hoppings_;
      std::vector<coeff_t> hopping_amplitudes_;
      std::vector<int> onsites_;
      std::vector<double> onsite_potentials_;
      std::vector<std::pair<int,int>> szszs_;
      std::vector<double> szsz_amplitudes_;
      std::vector<std::pair<int,int>> exchanges_;
      std::vector<coeff_t> exchange_amplitudes_;

      void initialize();

      state_t up_down_to_hole_table(const state_t& upspins, const state_t& downspin);
      state_t down_up_to_hole_table(const state_t& upspins, const state_t& downspin);

      inline int mpi_rank_of_spins(const state_t& spins) const
      { 
	state_t x = spins;
	x = ((x >> 16) ^ x) * 0x45d9f3b;
	x = ((x >> 16) ^ x) * 0x45d9f3b;
	x = (x >> 16) ^ x;
	return x % mpi_size_;
      }

      Spinhalf<state_t> hs_upspins_;
      Spinhalf<state_t> hs_downspins_;
      Spinhalf<state_t> hs_holes_in_ups_;      
      Spinhalf<state_t> hs_holes_in_downs_;
      LinTable<Spinhalf<state_t>, uint64> indexing_holes_in_ups_;
      LinTable<Spinhalf<state_t>, uint64> indexing_holes_in_downs_;

      std::vector<state_t> my_upspins_;
      std::unordered_map<state_t, uint64> my_upspins_offset_;
      std::vector<state_t> my_downspins_;
      std::unordered_map<state_t, uint64> my_downspins_offset_;

    
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

      std::vector<state_t> downspins_i_send_forward_;
      std::vector<state_t> upspins_i_send_forward_;
      std::vector<state_t> downspins_i_send_back_;
      std::vector<state_t> upspins_i_send_back_;

      std::vector<state_t> downspins_table_;
      std::vector<state_t> upspins_table_;

    

    };
  

  }
}


#endif
