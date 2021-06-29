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

#ifndef HYDRA_DYNAMICS_CONTINUEDFRACTION_
#define HYDRA_DYNAMICS_CONTINUEDFRACTION_

#include <vector>
#include <complex>

namespace hydra { namespace dynamics {

    std::complex<double> continued_fraction(const std::complex<double> z, 
					    const std::vector<double>& alphas,
					    const std::vector<double>& betas, 
					    int level=0);

  }
}
#endif
