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

#ifndef HYDRA_MODELS_HUBBARDMODEL_
#define HYDRA_MODELS_HUBBARDMODEL_

#include <utility>
#include <vector>

#include <hydra/bases/basis_electron.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <lila/matrix.h>

namespace hydra {

class HubbardModel {

public:

  HubbardModel(int n_sites, QN qn);
  lila::Matrix<coeff_t> matrix(BondList bondlist, Couplings couplings) const;

  void apply(BondList bondlist, Couplings couplings,
	     const lila::Vector<coeff_t> &in_vec,
	     lila::Vector<coeff_t> &out_vec) const;

  qn_t apply_fermion(const lila::Vector<coeff_t> &state_before,
                     lila::Vector<coeff_t> &state_after, std::string type,
                     int site) const;

  lila::Matrix<double> szMatrix(int siteIndex) const;
  lila::Matrix<double> sPlusMatrix(int siteIndex) const;
  lila::Matrix<double> sMinusMatrix(int siteIndex) const;

  qn_t QN() const { return qn_; }
  int n_sites() const { return n_sites_; }
  int64 dim() const { return dim_; }

private:
  int n_sites_;
  int64 dim_;
  qn_t qn_;

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
};

  lila::Matrix<coeff_t> matrix(BondList bondlist, Couplings couplings, 
			       HubbardModel const& model);
  void apply(BondList bondlist, Couplings couplings,
	     Hu

HubbardModel model, BondList bondlist, Couplings couplings)




} // namespace hydra

#endif
