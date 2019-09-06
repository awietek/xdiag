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

#ifndef HYDRA_MODELS_HUBBARDMODELDETAIL_
#define HYDRA_MODELS_HUBBARDMODELDETAIL_

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/typedefs.h>

namespace hydra { namespace models {


    namespace hubbardmodeldetail {
      using hydra::operators::BondList;
      using hydra::operators::Couplings;

      template <class coeff_t>
      void set_hubbard_terms(BondList bondlist, Couplings couplings,
			     std::vector<std::pair<int, int>>& hoppings,
			     std::vector<coeff_t>& hopping_amplitudes,
			     std::vector<std::pair<int, int>>& currents,
			     std::vector<coeff_t>& current_amplitudes,
			     std::vector<std::pair<int, int>>& interactions,
			     std::vector<double>& interaction_strengths,
			     std::vector<int>& onsites,
			     std::vector<double>& onsite_potentials, 
			     double& U);
      
      template <class coeff_t>
      void set_hoppings
	(BondList bondlist, Couplings couplings,
	 std::vector<std::pair<int, int>>& hoppings,
	 std::vector<coeff_t>& hopping_amplitudes);
      
      template <class coeff_t>
      void set_currents
      (BondList bondlist, Couplings couplings,
       std::vector<std::pair<int, int>>& currents,
       std::vector<coeff_t>& current_amplitudes);

      void set_interactions
      (BondList bondlist, Couplings couplings,
       std::vector<std::pair<int, int>>& interactions,
       std::vector<double>& interaction_strengths);
      
      void set_onsites
      (BondList bondlist, Couplings couplings,
       std::vector<int>& onsites,
       std::vector<double>& onsite_potentials);
      
      void set_U
      (Couplings couplings, double& U);
      
    }
  }
}

#endif
