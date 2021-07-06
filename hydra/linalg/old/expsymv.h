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

#ifndef LILA_SPARSE_EXPSYMV_H_
#define LILA_SPARSE_EXPSYMV_H_

#include "../matrix.h"
#include "../matrixfunction.h"
#include "../tmatrix.h"
#include "../vector.h"
#include "../common.h"

namespace lila {
    
  template <class multiply_f, class vector_t>
  real_t<typename vector_t::coeff_type>
  ExpSymVInplace(multiply_f A, vector_t& X, 
		 typename vector_t::coeff_type alpha, 
		 double precision, bool shift=false)
  {
    using coeff_type = typename vector_t::coeff_type;
    using real_type = real_t<coeff_type>;
    using tmatrix_t = Tmatrix<real_type>;
    
    auto norm = Norm(X);

    auto converged = 
      [precision, alpha, norm](const tmatrix_t& tmat, real_type beta) { 
      // return LanczosConvergedEigenvalues(tmat, beta, 1, (real_type)precision);
      return LanczosConvergedTimeEvolution(tmat, beta, alpha, (real_type)precision,
					   1000, norm);
    };
    auto v0 = X;
    auto res_first = Lanczos(A, v0, converged);
    v0.clear();
    v0.shrink_to_fit();

    auto tmat = Matrix<real_type>(res_first.tmatrix);
    int n_iterations = tmat.nrows();

    real_type e0 = EigenvaluesSym(tmat)(0);
    
    // Cast to complex (TODO: make this generic)
    auto tmatc = Zeros<coeff_type>(n_iterations, n_iterations);
    for (auto i : tmat.rows())
      for (auto j : tmat.cols())
	tmatc(i,j) = (coeff_type)tmat(i,j);

    if (shift)
      {
	for (auto i : tmat.rows())
	  tmatc(i,i) -= e0;
      }
    
    auto texp = ExpM(alpha * tmatc);
    std::vector<Vector<coeff_type>> linear_combinations = {texp.col(0)};

    auto converged_fixed = 
      [n_iterations](const tmatrix_t& tmat, real_type beta) {
      (void)beta;
      return LanczosConvergedFixed(tmat, n_iterations);
    };
    auto res = Lanczos(A, X, converged_fixed, linear_combinations);
    X = res.vectors[0] * (coeff_type)norm;
    return e0;
  }


  template <class multiply_f, class vector_t>
  vector_t ExpSymV(multiply_f A, const vector_t& X, 
		   typename vector_t::coeff_type alpha, double precision)
  {
    auto v0 = X;
    ExpSymVInplace(A, v0, alpha, precision);
    return v0;
  }

      
}

#endif
