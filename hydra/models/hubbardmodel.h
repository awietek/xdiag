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

#include <lila/matrix.h>
#include <hydra/hilbertspaces/hubbard.h>

namespace hydra { namespace models {
    
    /*!
      Class to generate representations of the Hubbard models
    */
    class HubbardModel {

    public:

      /*! 
	Defines a Hubbard model given number of sites and pairs of 
	neighboring sites with hoppings
      */
      HubbardModel(const int& n_sites, 
		   const std::vector<std::pair<int, int>> neighbors);
      
      ///  Return all possible quantum numbers
      std::vector<hilbertspaces::hubbard_qn> quantumnumbers();
      
      /*!
	returns a lila::Matrix of the Hubbard model given t, U, and
	the quantum numbers (n_upspins, n_downspins).

	Usage:
	@code
	#include <hydra/models/hubbardmodel.h>

	using hydra::hilbertspaces::hubbard_qn;
	using hydra::models::HubbardModel;

	int n_upspins = 4;
	int n_downspins = 4;
	int n_sites = 8;
	std::vector<std::pair<int, int>> hoppings;
	for (int i = 0; i < n_sites; ++i)
	  hoppings.push_back({i, (i + 1) % n_sites});
	hubbard_qn qn = {nup, ndown}; 
	auto model = HubbardModel(n_sites, hoppings);
	auto hamilton = model.matrix(t, U, qn);
	@endcode
      */
      lila::Matrix<double> matrix(const double& t, const double& U,
				  const hilbertspaces::hubbard_qn& qn) const;

      lila::Vector<double> apply_fermion(const lila::Vector<double>& state_before, 
					 hilbertspaces::hubbard_qn& qn, 
					 std::string type, int site) const;

    private:
      int n_sites_;
      std::vector<std::pair<int, int>> neighbors_;
    };
    
  }
}

#endif
