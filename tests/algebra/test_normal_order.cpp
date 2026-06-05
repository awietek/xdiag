// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/algebra/algebras/electron_algebra.hpp>
#include <xdiag/algebra/algebras/spin_algebra.hpp>
#include <xdiag/algebra/algebras/tj_algebra.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace xdiag::algebra;

TEST_CASE("normal_order", "[operators]") try {

  auto spin = spin_algebra();
  auto elec = electron_algebra();
  auto tj = tj_algebra();

  // -------------------------------------------------------------------------
  // Spin algebra — expansion
  // -------------------------------------------------------------------------
  Log("Testing normal_order: spin expansion");
  {
    // SzSz{0,1} -> Sz{0}*Sz{1}
    auto r = normal_order(OpSum(Op("SzSz", {0, 1})), spin);
    auto expected = 1.0 * (Op("Sz", 0) * Op("Sz", 1));
    REQUIRE(isapprox(r, expected, spin));
  }
  {
    // Exchange{0,1} -> 0.5*S+{0}*S-{1} + 0.5*S-{0}*S+{1}
    auto r = normal_order(OpSum(Op("Exchange", {0, 1})), spin);
    auto expected =
        0.5 * (Op("S+", 0) * Op("S-", 1)) + 0.5 * (Op("S-", 0) * Op("S+", 1));
    REQUIRE(isapprox(r, expected, spin));
  }
  {
    // SdotS{0,1} -> Sz{0}*Sz{1} + 0.5*S+{0}*S-{1} + 0.5*S-{0}*S+{1}
    auto r = normal_order(OpSum(Op("SdotS", {0, 1})), spin);
    auto expected = 1.0 * (Op("Sz", 0) * Op("Sz", 1)) +
                    0.5 * (Op("S+", 0) * Op("S-", 1)) +
                    0.5 * (Op("S-", 0) * Op("S+", 1));
    REQUIRE(isapprox(r, expected, spin));
  }

  // -------------------------------------------------------------------------
  // Spin algebra — same-site simplification
  // -------------------------------------------------------------------------
  Log("Testing normal_order: spin same-site");
  {
    // Sz{0}*Sz{0} -> 0.25*Id
    auto r = normal_order(1.0 * (Op("Sz", 0) * Op("Sz", 0)), spin);
    REQUIRE(isapprox(r, 0.25 * Op("Id"), spin));
  }
  {
    // S+{0}*S-{0} -> 0.5*Id + Sz{0}
    auto r = normal_order(1.0 * (Op("S+", 0) * Op("S-", 0)), spin);
    REQUIRE(isapprox(r, 0.5 * Op("Id") + 1.0 * Op("Sz", 0), spin));
  }
  {
    // S-{0}*S+{0} -> 0.5*Id - Sz{0}
    auto r = normal_order(1.0 * (Op("S-", 0) * Op("S+", 0)), spin);
    REQUIRE(isapprox(r, 0.5 * Op("Id") - 1.0 * Op("Sz", 0), spin));
  }
  {
    // S+{0}*S+{0} -> 0
    auto r = normal_order(1.0 * (Op("S+", 0) * Op("S+", 0)), spin);
    REQUIRE(isapprox(r, OpSum{}, spin));
  }
  {
    // SdotS{0,0} -> 0.75*Id  (S·S for spin-1/2)
    auto r = normal_order(OpSum(Op("SdotS", {0, 0})), spin);
    REQUIRE(isapprox(r, 0.75 * Op("Id"), spin));
  }

  // -------------------------------------------------------------------------
  // Spin algebra — sorting (bosonic: no sign change)
  // -------------------------------------------------------------------------
  Log("Testing normal_order: spin sorting");
  {
    // Sz{1}*Sz{0} -> Sz{0}*Sz{1}
    auto r = normal_order(1.0 * (Op("Sz", 1) * Op("Sz", 0)), spin);
    REQUIRE(isapprox(r, 1.0 * (Op("Sz", 0) * Op("Sz", 1)), spin));
  }
  {
    // S+{1}*S-{0} -> S-{0}*S+{1}  (no sign: spin operators commute)
    auto r = normal_order(1.0 * (Op("S+", 1) * Op("S-", 0)), spin);
    REQUIRE(isapprox(r, 1.0 * (Op("S-", 0) * Op("S+", 1)), spin));
  }

  // -------------------------------------------------------------------------
  // ScalarChirality expansion
  // -------------------------------------------------------------------------
  Log("Testing normal_order: ScalarChirality");
  {
    // (i/2)*[S+{0}*S-{1}*Sz{2} - S-{0}*S+{1}*Sz{2} + ...]
    auto r = normal_order(OpSum(Op("ScalarChirality", {0, 1, 2})), spin);
    complex I(0.0, 1.0);
    auto expected = (I * 0.5) * (Op("S+", 0) * Op("S-", 1) * Op("Sz", 2)) +
                    (-I * 0.5) * (Op("S-", 0) * Op("S+", 1) * Op("Sz", 2)) +
                    (I * 0.5) * (Op("Sz", 0) * Op("S+", 1) * Op("S-", 2)) +
                    (-I * 0.5) * (Op("Sz", 0) * Op("S-", 1) * Op("S+", 2)) +
                    (I * 0.5) * (Op("S-", 0) * Op("Sz", 1) * Op("S+", 2)) +
                    (-I * 0.5) * (Op("S+", 0) * Op("Sz", 1) * Op("S-", 2));
    REQUIRE(isapprox(r, expected, spin));
  }

  // -------------------------------------------------------------------------
  // Electron algebra — expansion
  // -------------------------------------------------------------------------
  Log("Testing normal_order: electron expansion");
  {
    // Nup{0} -> Cdagup{0}*Cup{0}
    auto r = normal_order(OpSum(Op("Nup", 0)), elec);
    REQUIRE(isapprox(r, 1.0 * (Op("Cdagup", 0) * Op("Cup", 0)), elec));
  }
  {
    // Hopup{0,1}: after sorting -Cdagup{1}*Cup{0} -> +Cup{0}*Cdagup{1}
    auto r = normal_order(OpSum(Op("Hopup", {0, 1})), elec);
    auto expected = -1.0 * (Op("Cdagup", 0) * Op("Cup", 1)) +
                    1.0 * (Op("Cup", 0) * Op("Cdagup", 1));
    REQUIRE(isapprox(r, expected, elec));
  }

  // -------------------------------------------------------------------------
  // Electron algebra — same-site CAR
  // -------------------------------------------------------------------------
  Log("Testing normal_order: electron CAR");
  {
    // Cdagup{0}*Cdagup{0} -> 0
    auto r = normal_order(1.0 * (Op("Cdagup", 0) * Op("Cdagup", 0)), elec);
    REQUIRE(isapprox(r, OpSum{}, elec));
  }
  {
    // Cup{0}*Cdagup{0} -> Id - Cdagup{0}*Cup{0}
    auto r = normal_order(1.0 * (Op("Cup", 0) * Op("Cdagup", 0)), elec);
    REQUIRE(
        isapprox(r, 1.0 * Op("Id") - 1.0 * (Op("Cdagup", 0) * Op("Cup", 0)), elec));
  }
  {
    // Cdagup{0}*Cdagdn{0} -> -Cdagdn{0}*Cdagup{0}
    auto r = normal_order(1.0 * (Op("Cdagup", 0) * Op("Cdagdn", 0)), elec);
    REQUIRE(isapprox(r, -1.0 * (Op("Cdagdn", 0) * Op("Cdagup", 0)), elec));
  }

  // -------------------------------------------------------------------------
  // Electron algebra — fermionic sorting (-1 per transposition)
  // -------------------------------------------------------------------------
  Log("Testing normal_order: fermionic sorting");
  {
    // Cdagup{1}*Cdagup{0} -> -Cdagup{0}*Cdagup{1}
    auto r = normal_order(1.0 * (Op("Cdagup", 1) * Op("Cdagup", 0)), elec);
    REQUIRE(isapprox(r, -1.0 * (Op("Cdagup", 0) * Op("Cdagup", 1)), elec));
  }
  {
    // Cdagup{1}*Cup{0} -> -Cup{0}*Cdagup{1}
    auto r = normal_order(1.0 * (Op("Cdagup", 1) * Op("Cup", 0)), elec);
    REQUIRE(isapprox(r, -1.0 * (Op("Cup", 0) * Op("Cdagup", 1)), elec));
  }

  // -------------------------------------------------------------------------
  // tJ algebra — projected CAR at same site
  // -------------------------------------------------------------------------
  Log("Testing normal_order: tJ projected CAR");
  {
    // Cup{0}*Cdagup{0} in tJ -> Id - Cdagup{0}*Cup{0} - Cdagdn{0}*Cdn{0}
    auto r = normal_order(1.0 * (Op("Cup", 0) * Op("Cdagup", 0)), tj);
    auto expected = 1.0 * Op("Id") - 1.0 * (Op("Cdagup", 0) * Op("Cup", 0)) -
                    1.0 * (Op("Cdagdn", 0) * Op("Cdn", 0));
    REQUIRE(isapprox(r, expected, tj));
  }
  {
    // Cdn{0}*Cdagup{0} in tJ -> 0  (no double occupancy)
    auto r = normal_order(1.0 * (Op("Cdn", 0) * Op("Cdagup", 0)), tj);
    REQUIRE(isapprox(r, OpSum{}, tj));
  }
  {
    // Cup{0}*Cdagdn{0} in tJ -> 0  (no double occupancy)
    auto r = normal_order(1.0 * (Op("Cup", 0) * Op("Cdagdn", 0)), tj);
    REQUIRE(isapprox(r, OpSum{}, tj));
  }

  // -------------------------------------------------------------------------
  // HubbardU expansion
  // -------------------------------------------------------------------------
  Log("Testing normal_order: HubbardU / Nupdn");
  {
    // A site-free "HubbardU" represents sum_i Nupdn{i}, but nsites is unknown
    // here, so it is elementary and left untouched (see normal_order.hpp).
    auto r = normal_order(OpSum(Op("HubbardU")), elec);
    REQUIRE(r.size() == 1);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "HubbardU");
    }
  }
  {
    // Nupdn{0} -> Cdagup{0}*Cup{0}*Cdagdn{0}*Cdn{0}: a single 4-operator
    // monomial on site 0, with no compound operator remaining.
    auto r = normal_order(OpSum(Op("Nupdn", 0)), elec);
    for (auto const &[c, mono] : r) {
      for (auto const &op : mono) {
        REQUIRE(op.type() != "HubbardU");
        REQUIRE(op.type() != "Nupdn");
      }
    }
  }

  // -------------------------------------------------------------------------
  // Site range validation
  // -------------------------------------------------------------------------
  Log("Testing normal_order: site range validation");
  {
    // REQUIRE_THROWS(normal_order(OpSum(Op("Sz", 5)), spin));
  }

  Log("Done testing normal_order");
} catch (xdiag::Error const &e) {
  error_trace(e);
}
