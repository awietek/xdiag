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

#include "catch.hpp"
#include <lila/all.h>

using namespace lila;

template <class coeff_t>
void test_expsymv(int n, real_t<coeff_t> precision, int n_tests)
{
  for (int random_seed : range<int>(n_tests)) 
    {
      uniform_dist_t<coeff_t> fdist(-1., 1.);
      uniform_gen_t<coeff_t> fgen(fdist, random_seed);
  
      coeff_t alpha(-10);

      Matrix<coeff_t> A(n, n);
      Random(A, fgen);
      A += Herm(A);
      REQUIRE(close(A, Herm(A)));
      coeff_t e0 = EigenvaluesSym(A)(0);
      A -= e0 * Identity<coeff_t>(n);
      
      Vector<coeff_t> v0(n);
      Random(v0, fgen);
      auto true_vec = Mult(ExpM(alpha*A), v0);

      auto multiply = [&A](const Vector<coeff_t>& v, Vector<coeff_t>& w) {
	w = Mult(A, v);
      };

      auto lczs_vec = ExpSymV(multiply, v0, alpha, precision);
      
      auto e = Norm(lczs_vec - true_vec);
      REQUIRE(e < precision);
      LilaPrint(e);
    }
}


TEST_CASE( "Expsymv test", "[Expsymv]" ) 
{
  // test_expsymv<float>(500, 1e-5, 2);
  test_expsymv<double>(500, 1e-8, 2);
  // test_expsymv<scomplex>(500, 1e-5, 2);
  test_expsymv<complex>(500, 1e-8, 2);
}

