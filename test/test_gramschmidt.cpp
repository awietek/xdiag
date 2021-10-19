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
void test_gramschmidt()
{
  using namespace lila;
  int n=20;
  int l=12;

  
  for (int seed : lila::range<int>(3)) {
  
    uniform_dist_t<coeff_t> fdist(-1., 1.);
    uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    std::vector<Vector<coeff_t>> vs;
    for (int i=0; i<l; ++i)
      {
	Vector<coeff_t> v(n);
	Random(v, fgen);
	vs.push_back(v);
      }

    auto es = gramschmidt(vs);
    for (int i=0; i<l; ++i)
    for (int j=0; j<l; ++j)
      {
	if (i==j) REQUIRE(close(lila::real(Dot(es[i], es[j])), (real_t<coeff_t>)1.));
	else REQUIRE(close(lila::real(Dot(es[i], es[j])), (real_t<coeff_t>)0.));
      }
    

  }
}

TEST_CASE( "gramschmidt test", "[gramschmidt]" ) {
  // test_gramschmidt<float>();
  test_gramschmidt<double>();
  // test_gramschmidt<std::complex<float>>();
  test_gramschmidt<std::complex<double>>();
}
