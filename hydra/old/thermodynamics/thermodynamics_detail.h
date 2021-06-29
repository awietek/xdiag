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

#ifndef HYDRA_THERMODYNAMICS_THERMODYNAMICS_DETAIL_
#define HYDRA_THERMODYNAMICS_THERMODYNAMICS_DETAIL_

#include <vector>
#include <cassert>


namespace hydra { namespace thermodynamics { namespace detail {

      void combine_thermodynamics
      (const std::vector<double> temperatures,
       const std::vector<double>& e0s,
       std::vector<std::vector<double>>& partitions_for_qn,
       std::vector<std::vector<double>>& energies_for_qn,
       std::vector<std::vector<double>>& quad_moments_for_qn,
       std::vector<double>& partitions,
       std::vector<double>& energies,
       std::vector<double>& quad_moments,
       std::vector<double>& specific_heats);

      double combined_quantity(const std::vector<double>& quantities, 
			       const std::vector<double>& e0s, double beta);


    }
  }
}

#endif
