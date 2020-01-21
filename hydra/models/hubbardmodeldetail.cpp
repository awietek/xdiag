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

#include <hydra/utils/complex.h>

#include "hubbardmodeldetail.h"

namespace hydra { namespace models {

    namespace hubbardmodeldetail {
      
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
			     std::vector<std::pair<int,int>> szszs,
			     std::vector<double> szsz_amplitudes,
			     std::vector<std::pair<int,int>> exchanges,
			     std::vector<coeff_t> exchange_amplitudes,
			     double& U)
      {
	set_hoppings(bondlist, couplings, hoppings, hopping_amplitudes);
	set_currents<coeff_t>(bondlist, couplings, currents, current_amplitudes);
	set_interactions(bondlist, couplings, interactions, 
			 interaction_strengths);
	set_onsites(bondlist, couplings, onsites, onsite_potentials);
	set_U(couplings, U);
	set_szszs(bondlist, couplings, szszs, szsz_amplitudes);
	set_exchanges<coeff_t>(bondlist, couplings, exchanges, exchange_amplitudes);

      }

      
      template <class coeff_t>
      void set_hoppings
	(BondList bondlist, Couplings couplings,
	 std::vector<std::pair<int, int>>& hoppings,
	 std::vector<coeff_t>& hopping_amplitudes)
      {
	BondList hopping_list = bondlist.bonds_of_type("HUBBARDHOP");
	for (auto bond : hopping_list)
	  {
	    int s1 = bond.sites()[0];
	    int s2 = bond.sites()[1];
	    hoppings.push_back({s1, s2});
	    coeff_t c = ForceReal<coeff_t>
	      (couplings[bond.coupling()], true, 
	       "Warning: deprecating imaginary part of hopping (real Hubbard)!");
	    hopping_amplitudes.push_back(c);
	  }
      }
      
      template <>
      void set_currents<double>
      (BondList bondlist, Couplings couplings,
       std::vector<std::pair<int, int>>& currents,
       std::vector<double>& current_amplitudes)
      {
	BondList current_list = bondlist.bonds_of_type("HUBBARDCURRENT");
	for (auto bond : current_list)
	{
	  int s1 = bond.sites()[0];
	  int s2 = bond.sites()[1];
	  currents.push_back({s1, s2});
	  assert( couplings.is_real(bond.coupling()) ); // current be real
	  double c = couplings.real(bond.coupling());
	  current_amplitudes.push_back(c);
	}
      }

      template <>
      void set_currents<complex>
      (BondList bondlist, Couplings couplings,
       std::vector<std::pair<int, int>>& currents,
       std::vector<complex>& current_amplitudes)
      {
	BondList current_list = bondlist.bonds_of_type("HUBBARDCURRENT");
	for (auto bond : current_list)
	{
	  int s1 = bond.sites()[0];
	  int s2 = bond.sites()[1];
	  currents.push_back({s1, s2});
	  assert( couplings.is_real(bond.coupling()) ); // current be real
	  complex c = complex(0., couplings.real(bond.coupling()));
	  current_amplitudes.push_back(c);
	}
      }

      void set_interactions
      (BondList bondlist, Couplings couplings,
       std::vector<std::pair<int, int>>& interactions,
       std::vector<double>& interaction_strengths)
      {
	BondList interaction_list = bondlist.bonds_of_type("HUBBARDV");
	for (auto bond : interaction_list)
	  {
	    int s1 = bond.sites()[0];
	    int s2 = bond.sites()[1];
	    interactions.push_back({s1, s2});
	    assert( couplings.is_real(bond.coupling()) );  // V be real
	    interaction_strengths.push_back(couplings.real(bond.coupling()));
	  }
      }
      
      void set_onsites
      (BondList bondlist, Couplings couplings,
       std::vector<int>& onsites,
       std::vector<double>& onsite_potentials)
      {
	BondList onsites_list = bondlist.bonds_of_type("HUBBARDMU");
	for (auto bond : onsites_list)
	  if (couplings.defined(bond.coupling()))
	    {
	      assert(bond.sites().size() == 1);
	      int s1 = bond.sites()[0];
	      onsites.push_back(s1);
	      assert( couplings.is_real(bond.coupling()) );  // mu be real
	      onsite_potentials.push_back(couplings.real(bond.coupling()));
	    }
      }
      
      void set_U
      (Couplings couplings, double& U)
      {
	if (couplings.defined("U"))
	  {
	    assert(couplings.is_real("U"));  // U be real
	    U = couplings.real("U");
	  }
	else U = 0.;
      }


      void set_szszs
	(BondList bondlist, Couplings couplings,
	 std::vector<std::pair<int, int>>& szszs,
	 std::vector<double>& szsz_amplitudes)
      {
	BondList szsz_list = bondlist.bonds_of_type("HEISENBERG");
	for (auto bond : szsz_list)
	  {
	    int s1 = bond.sites()[0];
	    int s2 = bond.sites()[1];
	    szszs.push_back({s1, s2});
	    double c = ForceReal<double>
	      (couplings[bond.coupling()], true, 
	       "Warning: deprecating imaginary part of Heisenberg (real Hubbard)!");
	    szsz_amplitudes.push_back(c);
	  }
      }

      template <class coeff_t>
      void set_exchanges
	(BondList bondlist, Couplings couplings,
	 std::vector<std::pair<int, int>>& exchanges,
	 std::vector<coeff_t>& exchange_amplitudes)
      {
	BondList exchange_list = bondlist.bonds_of_type("HEISENBERG");
	for (auto bond : exchange_list)
	  {
	    int s1 = bond.sites()[0];
	    int s2 = bond.sites()[1];
	    exchanges.push_back({s1, s2});
	    coeff_t c = ForceReal<coeff_t>
	      (couplings[bond.coupling()], true, 
	       "Warning: deprecating imaginary part of Heisenberg (real Hubbard)!");
	    exchange_amplitudes.push_back(c);
	  }
      }

      template void set_hubbard_terms<double>
      (BondList bondlist, Couplings couplings,
       std::vector<std::pair<int, int>>& hoppings,
       std::vector<double>& hopping_amplitudes,
       std::vector<std::pair<int, int>>& currents,
       std::vector<double>& current_amplitudes,
       std::vector<std::pair<int, int>>& interactions,
       std::vector<double>& interaction_strengths,
       std::vector<int>& onsites,
       std::vector<double>& onsite_potentials,
       std::vector<std::pair<int,int>> szszs,
       std::vector<double> szsz_amplitudes,
       std::vector<std::pair<int,int>> exchanges,
       std::vector<double> exchange_amplitudes,
       double& U);

      template void set_hubbard_terms<complex>
      (BondList bondlist, Couplings couplings,
       std::vector<std::pair<int, int>>& hoppings,
       std::vector<complex>& hopping_amplitudes,
       std::vector<std::pair<int, int>>& currents,
       std::vector<complex>& current_amplitudes,
       std::vector<std::pair<int, int>>& interactions,
       std::vector<double>& interaction_strengths,
       std::vector<int>& onsites,
       std::vector<double>& onsite_potentials,
       std::vector<std::pair<int,int>> szszs,
       std::vector<double> szsz_amplitudes,
       std::vector<std::pair<int,int>> exchanges,
       std::vector<complex> exchange_amplitudes,
       double& U);
  
    }
  }
}
