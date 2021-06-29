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

#ifndef HYDRA_THERMODYNAMICS_THERMODYNAMICS_EXACT_
#define HYDRA_THERMODYNAMICS_THERMODYNAMICS_EXACT_

#include <vector>
#include <string>
#include <algorithm>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra { namespace thermodynamics {
    using operators::BondList;
    using operators::Couplings;
    
    /*!
      Exact thermodynamics measurements using full diagonalization
    */
    struct thermodynamics_exact_result_t
    {
      std::vector<double> e0s;
      std::vector<lila::Vector<double>> eigenvalues;
      std::vector<std::vector<double>> partitions_for_qn;
      std::vector<std::vector<double>> energies_for_qn;
      std::vector<std::vector<double>> quad_moments_for_qn;
      std::vector<double> partitions;
      std::vector<double> energies;
      std::vector<double> quad_moments;
      std::vector<double> specific_heats;
    };

    template <class model_t>
    thermodynamics_exact_result_t
    thermodynamics_exact(const BondList& bondlist, 
			 const Couplings& couplings,
			 const std::vector<typename model_t::qn_t>& qns,
			 const std::vector<double>& temperatures,
			 bool keep_evs=false);
    
    template <class qn_t>
    void write_thermodynamics_exact
    (const std::vector<qn_t>& qns, const std::vector<double>& temperatures,
     const thermodynamics_exact_result_t& result, std::string outfile, 
     bool write_evs=false);

  }
}

#endif
