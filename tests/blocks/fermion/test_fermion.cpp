// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/matrices/matrix.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

using namespace xdiag;

TEST_CASE("fermion", "[fermion]") try {
  Log("Fermion anto-commutation relations");
  int nsites = 4;
  auto b = Fermion(nsites);

  for (int i=0; i<nsites; ++i){
    for (int j=0; j<nsites; ++j) {
      auto cdagi = matrix(Op("Cdag", i), b);
      auto cdagj = matrix(Op("Cdag", j), b);
      auto ci = matrix(Op("C", i), b);
      auto cj = matrix(Op("C", j), b);
      auto id = arma::eye(nsites, nsites);

      auto cdagcdag = cdagi * cdagj + cdagj * cdagi;
      auto cdagc = cdagi * cj + cj * cdagi;
      auto ccdag = ci * cdagj + cdagj * ci;
      auto cc = ci * cj + cj * ci;      
      
    }
  } 
} catch (xdiag::Error e) {
  error_trace(e);
}
