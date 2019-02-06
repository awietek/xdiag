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

#ifndef HYDRA_THERMODYNAMICS_THERMODYNAMICS_
#define HYDRA_THERMODYNAMICS_THERMODYNAMICS_

#include <vector>

#include <algorithm>

namespace hydra { namespace thermodynamics {

    /*!
      Combine a thermodynamic quantity from different blocks
      with different ground state energies at a given inverse
      temperature, according to the formula

      @f{equation}{
      A = \sum\limits_\rho e^{-\beta(e_{0,\rho} - e_0)} A_\rho
      @f}

      @param quantities the quantities in different sectors to combine,
      e.g. energy, partition function, quadratic moments
      @param e0s ground state energies of the different blocks
      @param beta inverse temperature
    */
    double combined_quantity(const std::vector<double>& quantities, 
			     const std::vector<double>& e0s, double beta);

  }
}
#endif
