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

template <class coeff_t>
void test_lobpcg()
{
  using namespace lila;

  int n=200;
  double precision = 10*rtol<coeff_t>::val();
  int max_iterations = 200;
  int num_eigenvalues = 5;

  // Test Lobpcg
  for (int random_seed : range<int>(2)) 
    {
  
      uniform_dist_t<coeff_t> dist(-1., 1.);
      uniform_gen_t<coeff_t> gen(dist, random_seed);
  
      Matrix<coeff_t> A(n, n);
      Random(A, gen);
      A += Herm(A);

      REQUIRE(close(A, Herm(A)));
      auto true_eigs = EigenvaluesSym(A);
  
      size_type n = A.nrows();
      std::vector<Vector<coeff_t>> vs;
      for (int i=0; i<num_eigenvalues; ++i)
	{
	  Vector<coeff_t> v(n);
	  Random(v, gen);
	  vs.push_back(v);
	}
   

      auto multiply = [&A](const Vector<coeff_t>& v, Vector<coeff_t>& w) {
	w = Mult(A, v);
      };
      auto res = Lobpcg(multiply, vs, precision, max_iterations);

      // LilaPrint(res.eigenvalues);
      // for (auto e : res.eigenvalues_history)
      //   printf("%.17g,", e(0) - true_eigs(0));
      // printf("\n");
      // for (auto e : res.eigenvalues_history)
      //   printf("%.17g,", e(1) - true_eigs(1));
      // printf("\n");
      // for (auto e : res.eigenvalues_history)
      //   printf("%.17g,", e(2) - true_eigs(2));
      // printf("\n");
      // for (auto e : res.eigenvalues_history)
      //   printf("%.17g,", e(3) - true_eigs(3));
      // printf("\n");
      // for (auto e : res.eigenvalues_history)
      //   printf("%.17g,", e(4) - true_eigs(4));
      // printf("\n");

      int idx=0;
      for(auto& v : res.eigenvectors)
	{
	  auto w = Mult(A, v);
	  coeff_t rayleigh = Dot(v,w) / Dot(v,v);
	  REQUIRE(std::abs(rayleigh - true_eigs(idx))<100*precision);
	  // LilaPrint(true_eigs(idx));
	  // LilaPrint(rayleigh);
	  ++idx;
	}	

    }
}

TEST_CASE( "Lobpcg test", "[Lobpcg]" ) {
  // test_lobpcg<float>();
  test_lobpcg<double>();
  // test_lobpcg<std::complex<float>>();
  test_lobpcg<std::complex<double>>();
}
