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

#ifndef HYDRA_MODELS_ELECTRON_
#define HYDRA_MODELS_ELECTRON_

#include <utility>
#include <vector>

#include <hydra/bases/basis_electron.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <lila/matrix.h>

namespace hydra {

  class Electron {
  public:
    Electron(int n_sites, QN qn);
    qn_t QN() const { return qn_; }
    int n_sites() const { return n_sites_; }
    int64 dim() const { return dim_; }
    
  private:
    int n_sites_;
    QN qn_;
    BasisElectron<uint32> basis_;
    IndexElectron<uint32, int64> index_;
    int64 dim_;
  };

  lila::Matrix<coeff_t> matrix(BondList bondlist, Couplings couplings, 
			       Electron const& model);
  
  void apply(BondList bondlist, Couplings couplings,
	     Electron const& block_in, lila::Vector<coeff_t> const& vec_in,
	     Electron const& block_out, lila::Vector<coeff_t> & vec_out)

} // namespace hydra

#endif
