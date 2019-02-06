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

#ifndef HYDRA_MODELS_HEISENBERGMODEL_
#define HYDRA_MODELS_HEISENBERGMODEL_

#include <utility>
#include <vector>

#include <hydra/symmetries/charactertable.h>

#include <lila/matrix.h>

namespace hydra { namespace models {

using symmetries::CharacterTable;

    /*!
      Class to generate representations of the Heisenberg models
    */
    class HeisenbergModel {
    public:
      /*! 
	Defines a Heisenberg model given number of sites and pairs of 
	neighboring sites with Heisenberg interactions
      */
      HeisenbergModel(int n_sites, 
		      const std::vector<std::pair<int, int>> neighbors);
      
      ///  Return all possible quantum numbers, i.e. number of upspins
      std::vector<int> quantumnumbers();
      
      /*!
	returns a lila::Matrix of the Heisenberg model given
	the quantum number (number of upspins).

	Usage:
	@code

	@endcode
      */
      lila::Matrix<double> matrix(double J, int qn) const;

      lila::Matrix<complex> matrix
      (double J, int qn, CharacterTable& character_table, 
       std::string representation_name) const;

    private:
      int n_sites_;
      std::vector<std::pair<int, int>> neighbors_;
    };
    
  }
}

#endif
