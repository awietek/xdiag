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
void test_bandlanczos()
{
  using namespace lila;

  int n=500;
  double precision = 1e-3;
  int max_iterations = 100;
  int n_lowest = 2;
  int num_eigenvalue = n_lowest;
  int n_bands = 3;

  // Test Bandlanczos
  for (int random_seed : range<int>(1)) {
  
    uniform_dist_t<coeff_t> fdist(-1., 1.);
    uniform_gen_t<coeff_t> fgen(fdist, random_seed);
  
    Matrix<coeff_t> A(n, n);
    Random(A, fgen);
    A += Herm(A);

    REQUIRE(close(A, Herm(A)));
    
    size_type dim = A.nrows();
    auto multiply = [&A](const Vector<coeff_t>& v, Vector<coeff_t>& w) {
      w = Mult(A, v);
    };
    auto lzs = BandLanczos<coeff_t, decltype(multiply)>
      (dim, random_seed, max_iterations, precision, num_eigenvalue, multiply, n_bands);
    auto lzs_eigenvalues = lzs.eigenvalues();
   
    auto true_eigs = EigenvaluesSym(A);
    // LilaPrint(true_eigs);
    // LilaPrint(lzs_eigenvalues.eigenvalues);
    // for (int i=0; i<lzs_eigenvalues.eigenvalues.size(); ++i)
    //   printf("%f %d\n", lzs_eigenvalues.eigenvalues(i), lzs_eigenvalues.multiplicity[i]);

    REQUIRE((int)lzs_eigenvalues.eigenvalues.size() == 
	    (int)lzs_eigenvalues.eigenvectors.size());
    REQUIRE((int)lzs_eigenvalues.multiplicity.size() == 
	    (int)lzs_eigenvalues.eigenvectors.size());    

    // for (int k = 0; k<=n_lowest; ++k)
    //   {
    // 	REQUIRE(std::abs(true_eigs(k) - lzs_eigenvalues.eigenvalues(k)) / 
    // 		std::abs(true_eigs(k)) < precision);
    //   }
  }
}



TEST_CASE( "Bandlanczos test", "[Bandlanczos]" ) {
  // test_bandlanczos<float>();
  test_bandlanczos<double>();
  // test_bandlanczos<std::complex<float>>();
  test_bandlanczos<std::complex<double>>();
}
