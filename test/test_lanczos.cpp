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
void test_lanczos(int n, int n_eigenvalue, 
		  real_t<coeff_t> precision, int n_tests,
		  bool verbose = false)
{
  for (int random_seed : range<int>(n_tests)) 
    {
      uniform_dist_t<coeff_t> fdist(-1., 1.);
      uniform_gen_t<coeff_t> fgen(fdist, random_seed);
  
      Matrix<coeff_t> A(n, n);
      Random(A, fgen);
      A += Herm(A);
      REQUIRE(close(A, Herm(A)));
      Vector<coeff_t> v0(n);
      Random(v0, fgen);
      auto multiply = [&A](const Vector<coeff_t>& v, Vector<coeff_t>& w) {
	w = Mult(A, v);
      };
    
      auto true_eigs = EigenvaluesSym(A);

      // Lanczos eigenvalue computation with precision on eigenvalues
      auto tmatrix = LanczosEigenvalues(multiply, v0, precision, 
					n_eigenvalue, "Eigenvalues").tmatrix;
      auto lczs_eigs = Eigenvalues(tmatrix);
      
      for (int i=0; i<=n_eigenvalue; ++i)
	{
	  real_t<coeff_t> e = std::abs(true_eigs(i) - lczs_eigs(i)) 
	    / std::abs(true_eigs(i));
	  if (verbose) printf("Eigenvalue %d: %g\n", i, e);
	  if (i==0) REQUIRE(e < precision * (real_t<coeff_t>)10.);
	}

      // Lanczos eigenvalue computation with precision on Ritz value
      auto tmatrix_ritz = LanczosEigenvalues(multiply, v0, precision, 
					     n_eigenvalue, "Ritz").tmatrix;
      auto lczs_eigs_ritz = Eigenvalues(tmatrix_ritz);
      for (int i=0; i<=n_eigenvalue; ++i)
	{
	  real_t<coeff_t> e = std::abs(true_eigs(i) - lczs_eigs_ritz(i)) 
	    / std::abs(true_eigs(i));
	  if (verbose) 	  printf("Ritz %d: %g\n", i, e);
	  if (i==0) CHECK(e <= precision * 10);
	}

      auto res = LanczosEigenvectors<Vector<coeff_t>>(multiply, v0, fgen, true, 
						      precision, {0}, "Ritz");
      auto& v = res.vectors[0];
      auto Hv = v;
      multiply(v, Hv);
      auto e = std::abs(Dot(v, Hv) - true_eigs(0));
      REQUIRE(e < precision * (real_t<coeff_t>)100.);
    }
}



TEST_CASE( "Lanczos test", "[Lanczos]" ) {
  bool verbose = false;
  // test_lanczos<float>(100, 0, 1e-5, 2, verbose);
  test_lanczos<double>(100, 0, 1e-12, 2, verbose);
  // test_lanczos<lila::scomplex>(100, 0, 1e-5, 2, verbose);
  test_lanczos<lila::complex>(100, 0, 1e-12, 2, verbose);

  // test_lanczos<float>(200, 1, 1e-5, 2, verbose);
  test_lanczos<double>(200, 1, 1e-12, 5, verbose);
  // test_lanczos<lila::scomplex>(200, 1, 1e-5, 2, verbose);
  test_lanczos<lila::complex>(200, 1, 1e-12, 2, verbose);

  // test_lanczos<float>(400, 2, 1e-5, 2, verbose);
  test_lanczos<double>(400, 2, 1e-12, 2, verbose);
  // test_lanczos<lila::scomplex>(400, 2, 1e-5, 2, verbose);
  test_lanczos<lila::complex>(400, 2, 1e-12, 2, verbose);
}
