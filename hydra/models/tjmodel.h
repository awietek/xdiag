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

#include <lila/matrix.h>
#include <hydra/hilbertspaces/hubbard.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra { namespace models {
    using operators::BondList;
    using operators::Couplings;

    /*!
      Class to generate representations of the TJ models
    */
    template <class coeff_t>
    class TJModel {

    public:
      using qn_t = hilbertspaces::hubbard_qn;

      /*! 
	Defines a TJ model given number of sites and pairs of 
	neighboring sites with hoppings
      */
      TJModel(BondList bondlist, Couplings couplings, qn_t qn);
      
      /*!
	returns a lila::Matrix of the TJ model given t, U, and
	the quantum numbers (n_upspins, n_downspins).
      */
      lila::Matrix<coeff_t> matrix() const;

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
      std::vector<std::pair<int,int>> szszs_;
      std::vector<double> szsz_amplitudes_;
      std::vector<std::pair<int,int>> exchanges_;
      std::vector<coeff_t> exchange_amplitudes_;
      double U_;
    };
    
  }
}

#endif
