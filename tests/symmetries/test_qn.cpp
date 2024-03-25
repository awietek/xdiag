#include "../catch.hpp"

#include <iostream>

#include "../blocks/electron/testcases_electron.h"
#include <hydra/symmetries/qn.h>


TEST_CASE("QN", "[symmetries]") {
  using namespace hydra;

  U1 g;
  QNum q1(g, 2);
  QNum q2(g, 3);
  QNum q3(g, 5);
  auto qn = q1 * q2;
  // HydraPrint(q1);
  // HydraPrint(q2);
  // HydraPrint(q3);
  // HydraPrint(qn);
  REQUIRE(qn == q3);

  int n_sites = 5;
  auto [space_group, irreps] =
      hydra::testcases::electron::get_cyclic_group_irreps(n_sites);
  QNum qq1(space_group, irreps[1]);
  QNum qq2(space_group, irreps[2]);
  QNum qq3(space_group, irreps[3]);
  auto qqn = qq1 * qq2;
  // HydraPrint(qq1);
  // HydraPrint(qq2);
  // HydraPrint(qq3);
  // HydraPrint(qqn);
  REQUIRE(qqn == qq3);
}
