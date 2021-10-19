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
#include <lila/allmpi.h>
#include <unordered_map>
#include <vector>

#include <hydra/states/state_electron.h>
#include <hydra/bases/basis_spinhalf.h>
#include <hydra/indexing/index_table.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class coeff_t, class bit_t = std_bit_t, class idx_t = std_idx_t>
class HubbardModelMPI {
public:
  using state_t = state_electron<bit_t>;

  HubbardModelMPI(BondList bondlist, Couplings couplings, qn_electron qn);

  void apply_hamiltonian(const lila::VectorMPI<coeff_t> &in_vec,
                         lila::VectorMPI<coeff_t> &out_vec,
                         bool verbose = false);

  qn_electron apply_fermion(const lila::VectorMPI<coeff_t> &in_vec,
                            lila::VectorMPI<coeff_t> &out_vec, std::string type,
                            int site);

  qn_electron qn() const { return qn_; }
  void set_qn(qn_electron qn);
  int n_sites() const { return n_sites_; }
  idx_t local_dim() const { return local_dim_; }
  idx_t dim() const { return dim_; }
  idx_t index(bit_t upspins, bit_t downspins) {
    return mpi_rank_of_spins(upspins) == mpi_rank_
               ? my_upspins_offset_[upspins] +
      index_dn_.index({downspins})
               : -1;
  }

  std::vector<bit_t> my_upspins() { return my_upspins_; }
  std::unordered_map<bit_t, idx_t> my_upspins_offset() {
    return my_upspins_offset_;
  }

private:
  int n_sites_;
  qn_electron qn_;

  std::vector<std::pair<int, int>> hoppings_;
  std::vector<coeff_t> hopping_amplitudes_;
  std::vector<std::pair<int, int>> currents_;
  std::vector<coeff_t> current_amplitudes_;
  std::vector<std::pair<int, int>> interactions_;
  std::vector<double> interaction_strengths_;
  std::vector<int> onsites_;
  std::vector<double> onsite_potentials_;
  std::vector<std::pair<int, int>> szszs_;
  std::vector<double> szsz_amplitudes_;
  std::vector<std::pair<int, int>> exchanges_;
  std::vector<coeff_t> exchange_amplitudes_;

  double U_;

  // members below are all set by initialize
  void initialize();
  BasisSpinHalf<bit_t> basis_up_;
  BasisSpinHalf<bit_t> basis_dn_;
  IndexTable<BasisSpinHalf<bit_t>, idx_t> index_up_;
  IndexTable<BasisSpinHalf<bit_t>, idx_t> index_dn_;

  int mpi_rank_;
  int mpi_size_;

  idx_t dim_;
  idx_t local_dim_;
  idx_t max_local_dim_;
  idx_t min_local_dim_;
  idx_t local_dim_downspins_;

  std::vector<bit_t> my_upspins_;
  std::unordered_map<bit_t, idx_t> my_upspins_offset_;
  std::vector<bit_t> my_downspins_;
  std::unordered_map<bit_t, idx_t> my_downspins_offset_;

  // Determine process of a spin configuration
  inline int mpi_rank_of_spins(const bit_t &spins) const {
    // return (int)(std::hash<bit_t>{}(spins) % mpi_size_);
    bit_t x = spins;
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

  std::vector<bit_t> downspins_i_recv_forward_;
  std::vector<bit_t> upspins_i_recv_forward_;
  std::vector<int> downspins_i_send_forward_offsets_;
  std::vector<int> downspins_i_recv_forward_offsets_;

  idx_t sum_n_downspins_i_send_forward_;
  idx_t sum_n_downspins_i_recv_forward_;

  std::vector<int> n_upspins_i_send_back_;
  std::vector<int> n_upspins_i_recv_back_;
  std::vector<int> n_upspins_i_send_back_offsets_;
  std::vector<int> n_upspins_i_recv_back_offsets_;
  std::vector<bit_t> downspins_i_recv_back_;
  std::vector<bit_t> upspins_i_recv_back_;
  std::vector<int> upspins_i_send_back_offsets_;
  std::vector<int> upspins_i_recv_back_offsets_;

  idx_t sum_n_upspins_i_send_back_;
  idx_t sum_n_upspins_i_recv_back_;

  idx_t buffer_size_;
  std::vector<coeff_t> send_buffer_;
  std::vector<coeff_t> recv_buffer_;

  // std::vector<bit_t> downspins_i_send_forward_;
  // std::vector<bit_t> upspins_i_send_forward_;
  // std::vector<bit_t> downspins_i_send_back_;
  // std::vector<bit_t> upspins_i_send_back_;
};

} // namespace hydra

#endif
