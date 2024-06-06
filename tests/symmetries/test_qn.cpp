#include "../catch.hpp"

#include <iostream>

#include "../blocks/electron/testcases_electron.hpp"
#include <xdiag/symmetries/qn.hpp>


TEST_CASE("QN", "[symmetries]") {
  using namespace xdiag;

  U1 g;
  QNum q1(g, 2);
  QNum q2(g, 3);
  QNum q3(g, 5);
  auto qn = q1 * q2;
  // XDIAG_SHOW(q1);
  // XDIAG_SHOW(q2);
  // XDIAG_SHOW(q3);
  // XDIAG_SHOW(qn);
  REQUIRE(qn == q3);

  int n_sites = 5;
  auto [space_group, irreps] =
      xdiag::testcases::electron::get_cyclic_group_irreps(n_sites);
  QNum qq1(space_group, irreps[1]);
  QNum qq2(space_group, irreps[2]);
  QNum qq3(space_group, irreps[3]);
  auto qqn = qq1 * qq2;
  // XDIAG_SHOW(qq1);
  // XDIAG_SHOW(qq2);
  // XDIAG_SHOW(qq3);
  // XDIAG_SHOW(qqn);
  REQUIRE(qqn == qq3);
}
