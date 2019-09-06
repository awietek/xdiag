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

#ifndef HYDRA_THERMODYNAMICS_THERMODYNAMICS_TPQ_
#define HYDRA_THERMODYNAMICS_THERMODYNAMICS_TPQ_

#include <lila/all.h>

namespace hydra { namespace thermodynamics {
      
    struct thermodynamics_tpq_result_t
    {
      double e0;
      lila::Tmatrix<double> tmatrix;
      lila::Vector<double> eigenvalues;
      lila::Vector<double> temperatures;
      lila::Vector<double> partitions;
      lila::Vector<double> energies;
      lila::Vector<double> quad_moments;
    };

    template <class model_t, class vector_t, class logger_t=lila::Logger>
    thermodynamics_tpq_result_t
    thermodynamics_tpq(model_t& model, vector_t& v0,
		       std::vector<double> temperatures,
		       double mintemperature = 0.01, 
		       double precision = 1e-12,
		       int iters = 1000, 
		       logger_t& logger=lila::logger);
    
    void write_thermodynamics_tpq_outfile
    (const thermodynamics_tpq_result_t& res, std::ofstream& of, 
     std::string comment="");

    void write_thermodynamics_tpq_tmatrix
    (const thermodynamics_tpq_result_t& res, std::ofstream& of, 
     std::string comment="");

    
  }
}

#endif
