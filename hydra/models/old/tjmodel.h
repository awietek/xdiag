
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

#ifndef HYDRA_MODELS_TJMODEL_
#define HYDRA_MODELS_TJMODEL_

#include <utility>
#include <vector>

#include <hydra/qns/qn_tj.h>
#include <hydra/bases/basis_tj.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <lila/matrix.h>

namespace hydra {

template <class coeff_t_, class bit_t_ = std_bit_t, class idx_t_ = std_idx_t>
class TJModel {

public:
  using coeff_t = coeff_t_;
  using bit_t = bit_t_;
  using idx_t = idx_t_;  
  using state_t = state_tj<bit_t>;
  using qn_t = qn_tj;
  using basis_t = BasisTJ<bit_t>;
  using vector_t = lila::Vector<coeff_t>;
  
  TJModel(BondList bondlist, Couplings couplings, qn_t qn);

  lila::Matrix<coeff_t> matrix(bool ninj_term = false) const;

  lila::Matrix<coeff_t> szMatrix(int siteIndex) const;
  lila::Matrix<double> sPlusMatrix(int siteIndex) const;
  lila::Matrix<double> sMinusMatrix(int siteIndex) const;

  qn_t qn() const { return qn_; }
  void set_qn(qn_t qn);
  int n_sites() const { return n_sites_; }
  int64 dim() const { return dim_; }

private:
  int n_sites_;
  int64 dim_;
  qn_t qn_;

  std::vector<std::pair<int, int>> hoppings_;
  std::vector<coeff_t> hopping_amplitudes_;
  std::vector<int> onsites_;
  std::vector<double> onsite_potentials_;
  std::vector<std::pair<int, int>> szszs_;
  std::vector<double> szsz_amplitudes_;
  std::vector<std::pair<int, int>> exchanges_;
  std::vector<coeff_t> exchange_amplitudes_;
  double U_;
};

} // namespace hydra

#endif
