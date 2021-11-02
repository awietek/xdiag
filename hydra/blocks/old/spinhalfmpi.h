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

#ifndef HYDRA_MODELS_SPINHALFMPI_
#define HYDRA_MODELS_SPINHALFMPI_

#include <lila/allmpi.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/indexing/lintable.h>

namespace hydra { namespace models {

    using operators::BondList;
    using operators::Couplings;
    using indexing::LinTable;


    class SpinhalfMPI {
    public:
      SpinhalfMPI();

      template <class coeff_t>
      quantumnumber multiply(quantumnumber qn,
			     BondList const& bondlist,
			     Couplings const& couplings,
			     lila::VectorMPI<coeff_t> const& x,
			     lila::VectorMPI<coeff_t> & y,
			     coeff_t alpha = (coeff_t)1.0,
			     coeff_t beta  = (coeff_t)0.0);

      uint64 dim(int n_sites, quantumnumber qn) const;
      bool block_defined(int n_sites, quantumnumber qn) const;
      void create_block(int n_sites, quantumnumber qn);

    private:
      int mpi_rank_;
      int mpi_size_;


      multiplyIsing(double J, std::vector<std::pair<int,int> sites,
		    lila::VectorMPI<coeff_t> const& x,
		    lila::VectorMPI<coeff_t> & y,
				coeff_t alpha = (coeff_t)1.0,
				coeff_t beta  = (coeff_t)0.0)
      
      std::vector<int> n_sites_;
      std::map<int, std::vector<LinTable>> lintables_;

      std::vector<std::pair<int, int>> n_sites_n_upspins_;
      std::vector<std::pair<int, int>, uint64> dim_;
      std::vector<std::pair<int, int>, uint64> dim_local_;
      std::vector<std::pair<int, int>, int> n_prefix_;
      std::vector<std::pair<int, int>, int> n_postfix_;
      
      std::map<std::pair<int, int>, std::vector<uint64>> my_prefixes_;
      std::map<std::pair<int, int>, std::vector<uint64>> my_postfixes_;
      std::map<std::pair<int, int>, std::vector<uint64>> my_prefix_size_;
      std::map<std::pair<int, int>, std::vector<uint64>> my_postfix_size_;
      std::map<std::pair<int, int>, std::vector<uint64>> my_prefix_offset_;
      std::map<std::pair<int, int>, std::vector<uint64>> my_postfix_offset_;

    };
    
  }
}

#endif
