// Copyright 2018 Alexander Wietek - All Rights Reserved.
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

#ifndef LILA_ALGORITHM_GRAMSCHMIDT_H_
#define LILA_ALGORITHM_GRAMSCHMIDT_H_

#include <vector>

namespace lila {
  template <class vector_t>
  std::vector<vector_t> gramschmidt(std::vector<vector_t>& ws, 
				    bool iterate=true)
  {
    auto es = ws;

    // Iterations for higher precision on orthogonality
    int n_iters = (iterate) ? 3 : 1;
    for (int iter=0; iter<n_iters; ++iter)
      {
	
	// Modified Gram-Schmidt
	for (int i=0; i<(int)ws.size(); ++i)
	  {
	    for (int k=0; k<i; ++k)
	      {
		auto rki = Dot(es[k], es[i]);
		es[i] -= rki * es[k];
	      }
	    auto rii = (typename vector_t::coeff_type)Norm(es[i]);
	    es[i] /= rii;
	  }

      }
    return es;   
  }

  template <class vector_t>
  bool has_full_rank(std::vector<vector_t>& ws)
  {
    using coeff_t = typename vector_t::coeff_type;

    int size = (int)ws.size();
    auto ovlps = lila::Zeros<coeff_t>(size, size);
    for (int i=0; i<size; ++i)
      {
    	ovlps(i, i) = lila::Dot(ws[i], ws[i]);
    	for (int j=i+1; j<size; ++j)
    	  {
    	    ovlps(i, j) = lila::Dot(ws[i], ws[j]);
    	    ovlps(j, i) = lila::conj(ovlps(i, j));
    	  }
      }
    auto eigs = EigenvaluesSym(ovlps);
    for (auto e : eigs)
      if (close<coeff_t>(e, (coeff_t)0.)) return false;
    return true;
  }


}

#endif
