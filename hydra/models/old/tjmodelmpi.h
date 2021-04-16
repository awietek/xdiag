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

#include <unordered_map>
#include <vector>

#include <hydra/qns/qn_tj.h>

#include <hydra/states/state_tj.h>

#include <hydra/bases/basis_spinhalf.h>
#include <hydra/bases/basis_tj.h>

#include <hydra/indexing/index_table.h>
#include <hydra/indexing/lintable.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <lila/allmpi.h>

namespace hydra {

template <class coeff_t_, class bit_t_ = std_bit_t, class idx_t_ = std_idx_t>
class TJModelMPI {
public:
  using coeff_t = coeff_t_;
  using bit_t = bit_t_;
  using idx_t = idx_t_;
  using state_t = state_tj<bit_t>;
  using qn_t = qn_tj;
  using basis_t = BasisTJ<bit_t>;
  using vector_t = lila::VectorMPI<coeff_t>;

  TJModelMPI(BondList bondlist, Couplings couplings, qn_tj qn);

  void apply_hamiltonian(lila::VectorMPI<coeff_t> const &in_vec,
                         lila::VectorMPI<coeff_t> &out_vec,
                         bool verbose = false);

  void apply_sz(const lila::VectorMPI<coeff_t> &in_vec,
                lila::VectorMPI<coeff_t> &out_vec, int site);

  qn_tj qn() const { return qn_; }
  int n_sites() const { return n_sites_; }
  idx_t local_dim() const { return local_dim_; }
  idx_t dim() const { return dim_; }

  template <class functor_t> void transverse(functor_t functor) const;

private:
  const int n_sites_;
  qn_tj qn_;

  int mpi_rank_;
  int mpi_size_;

  // Dimensions
  idx_t dim_;
  idx_t local_dim_;
  idx_t max_local_dim_;
  idx_t min_local_dim_;
  idx_t local_dim_downspins_;

  std::vector<std::pair<int, int>> hoppings_;
  std::vector<coeff_t> hopping_amplitudes_;
  std::vector<int> onsites_;
  std::vector<double> onsite_potentials_;
  std::vector<std::pair<int, int>> szszs_;
  std::vector<double> szsz_amplitudes_;
  std::vector<std::pair<int, int>> exchanges_;
  std::vector<coeff_t> exchange_amplitudes_;

  void initialize();

  bit_t up_down_to_hole_table(const bit_t &upspins, const bit_t &downspin);
  bit_t down_up_to_hole_table(const bit_t &upspins, const bit_t &downspin);

  inline int mpi_rank_of_spins(const bit_t &spins) const {
    bit_t x = spins;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x % mpi_size_;
  }

  BasisSpinHalf<bit_t> hs_upspins_;
  BasisSpinHalf<bit_t> hs_downspins_;
  BasisSpinHalf<bit_t> hs_holes_in_ups_;
  BasisSpinHalf<bit_t> hs_holes_in_downs_;
  LinTable<BasisSpinHalf<bit_t>, idx_t> indexing_holes_in_ups_;
  LinTable<BasisSpinHalf<bit_t>, idx_t> indexing_holes_in_downs_;

  std::vector<bit_t> my_upspins_;
  std::unordered_map<bit_t, idx_t> my_upspins_offset_;
  std::vector<bit_t> my_downspins_;
  std::unordered_map<bit_t, idx_t> my_downspins_offset_;

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

  std::vector<bit_t> downspins_i_send_forward_;
  std::vector<bit_t> upspins_i_send_forward_;
  std::vector<bit_t> downspins_i_send_back_;
  std::vector<bit_t> upspins_i_send_back_;

  std::vector<bit_t> downspins_table_;
  std::vector<bit_t> upspins_table_;
};

template <class coeff_t, class bit_t, class idx_t>
template <class functor_t>
void TJModelMPI<coeff_t, bit_t, idx_t>::transverse(functor_t functor) const {
  idx_t n_hole_configurations = hs_holes_in_ups_.size();
  for (bit_t const &upspins : my_upspins_) {
    idx_t upspin_offset = my_upspins_offset_.at(upspins);
    for (idx_t idx = upspin_offset; idx < upspin_offset + n_hole_configurations;
         ++idx) {
      bit_t downspins = downspins_table_[idx];
      functor({upspins, downspins}, idx);
    }
  }
}

} // namespace hydra

#endif
