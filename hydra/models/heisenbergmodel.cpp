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

#include <lila/special.h>

#include <hydra/symmetries/charactertable.h>
#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/indexing/indexspinhalf.h>
#include <hydra/indexing/indexsymmetrized.h>
#include <hydra/utils/range.h>
#include <hydra/utils/bitops.h>

#include "heisenbergmodel.h"

namespace hydra { namespace models {
    
    HeisenbergModel::HeisenbergModel
    (int n_sites, const std::vector<std::pair<int, int>> neighbors)
      : n_sites_(n_sites),
	neighbors_(neighbors)
    {}

    std::vector<int> HeisenbergModel::quantumnumbers()
    {
      std::vector<int> qns;
      for (int n_upspins = 0; n_upspins <= n_sites_; ++n_upspins)
	qns.push_back(n_upspins);
      return qns;
    }

    lila::Matrix<double> HeisenbergModel::matrix(double J, int qn) const
    {
      using hilbertspaces::Spinhalf;
      using indexing::IndexSpinhalf;
      using utils::range;
      using utils::gbit;
      using utils::popcnt;
	
      Spinhalf<uint64> hs(n_sites_, qn);
      IndexSpinhalf<uint64, uint64> indexing(hs);
      int dim = indexing.size();
      lila::Matrix<double> hamilton(dim, dim);
      lila::Zeros(hamilton);

      for (auto pair : neighbors_)
	{
	  int s1 = pair.first; 
	  int s2 = pair.second;
	  
	  // Apply Heisenberg operator on sites s1, s2
	  uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);
	  for (int idx : range<>(indexing.size()))
	    {
	      auto state = indexing.state(idx);
		
	      if (gbit(state, s1) == gbit(state, s2))
		hamilton(idx, idx) += J / 4.;  // Ising
	      else
		{
		  hamilton(idx, idx) -= J / 4.; // Ising
		  
		  //Exchange term
		  auto new_state = state ^ flipmask;
		  int new_idx = indexing.index(new_state);
		  hamilton(new_idx, idx) += J / 2.;
		}
	    }
	}
      return hamilton;
    }

    lila::Matrix<complex> HeisenbergModel::matrix
    (double J, int qn, CharacterTable& character_table, 
     std::string representation_name) const
    {
      using hilbertspaces::Spinhalf;
      using indexing::IndexSymmetrized;
      using utils::range;
      using utils::gbit;
      using utils::popcnt;
	
      Spinhalf<uint64> hs(n_sites_, qn);
      IndexSymmetrized<Spinhalf<uint64>> indexing
	(hs, character_table, representation_name);
      std::vector<complex> characters = 
	character_table.characters(representation_name);
      int dim = indexing.size();
      lila::Matrix<complex> hamilton(dim, dim);
      lila::Zeros(hamilton);

      for (auto pair : neighbors_)
	{
	  int s1 = pair.first; 
	  int s2 = pair.second;
	  
	  // Apply Heisenberg operator on sites s1, s2
	  uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);
	  for (int idx : range<>(indexing.size()))
	    {
	      auto state = indexing.state(idx);

	      if (gbit(state, s1) == gbit(state, s2))
		hamilton(idx, idx) += J / 4.;  // Ising
	      else
	      	{
	      	  hamilton(idx, idx) -= J / 4.; // Ising

	      	  //Exchange term
	      	  auto new_state = state ^ flipmask;
	      	  auto representative = new_state; 
	      	  int n_sym = indexing.find_representative(&representative);
	      	  int new_idx = indexing.index(representative);

	      	  if (new_idx != -1)
	      	    {
	      	      complex character = characters[n_sym];
	      	      complex coeff = indexing.norm(new_idx) / indexing.norm(idx) 
	      	      	* character * J / 2.;
	      	      hamilton(new_idx, idx) += coeff;
	      	    }	
	      	}
	    }
	}

      return hamilton;
    }
    
  }
}
