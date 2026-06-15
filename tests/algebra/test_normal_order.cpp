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

  // -------------------------------------------------------------------------
  // Spin algebra — expansion
  // -------------------------------------------------------------------------
  Algebra alg = spin_algebra(2);
  Log("Testing normal_order: spin expansion");
  {
    // SzSz{0,1} -> Sz{0}*Sz{1}
    auto r = normal_order(OpSum(Op("SzSz", {0, 1})), alg);
    auto expected = 1.0 * (Op("Sz", 0) * Op("Sz", 1));
    REQUIRE(isapprox(r, expected, alg));
  }
  {
    // Exchange{0,1} -> 0.5*S+{0}*S-{1} + 0.5*S-{0}*S+{1}
    auto r = normal_order(OpSum(Op("Exchange", {0, 1})), alg);
    auto expected =
        0.5 * (Op("S+", 0) * Op("S-", 1)) + 0.5 * (Op("S-", 0) * Op("S+", 1));
    REQUIRE(isapprox(r, expected, alg));
  }
  {
    // SdotS{0,1} -> Sz{0}*Sz{1} + 0.5*S+{0}*S-{1} + 0.5*S-{0}*S+{1}
    auto r = normal_order(OpSum(Op("SdotS", {0, 1})), alg);
    auto expected = 1.0 * (Op("Sz", 0) * Op("Sz", 1)) +
                    0.5 * (Op("S+", 0) * Op("S-", 1)) +
                    0.5 * (Op("S-", 0) * Op("S+", 1));
    REQUIRE(isapprox(r, expected, alg));
  }

  // -------------------------------------------------------------------------
  // Spin algebra — same-site simplification
  // -------------------------------------------------------------------------
  Log("Testing normal_order: spin same-site");
  alg = spin_algebra(1);
  {
    // Sz{0}*Sz{0} -> 0.25*Id
    auto r = normal_order(1.0 * (Op("Sz", 0) * Op("Sz", 0)), alg);
    REQUIRE(isapprox(r, 0.25 * Op("Id"), alg));
  }
  {
    // S+{0}*S-{0} -> 0.5*Id + Sz{0}
    auto r = normal_order(1.0 * (Op("S+", 0) * Op("S-", 0)), alg);
    REQUIRE(isapprox(r, 0.5 * Op("Id") + 1.0 * Op("Sz", 0), alg));
  }
  {
    // S-{0}*S+{0} -> 0.5*Id - Sz{0}
    auto r = normal_order(1.0 * (Op("S-", 0) * Op("S+", 0)), alg);
    REQUIRE(isapprox(r, 0.5 * Op("Id") - 1.0 * Op("Sz", 0), alg));
  }
  {
    // S+{0}*S+{0} -> 0
    auto r = normal_order(1.0 * (Op("S+", 0) * Op("S+", 0)), alg);
    REQUIRE(isapprox(r, OpSum{}, alg));
  }
  {
    // SdotS{0,0} -> 0.75*Id  (S·S for spin-1/2)
    auto r = normal_order(OpSum(Op("SdotS", {0, 0})), alg);
    REQUIRE(isapprox(r, 0.75 * Op("Id"), alg));
  }

  // -------------------------------------------------------------------------
  // Spin algebra — sorting (bosonic: no sign change)
  // -------------------------------------------------------------------------
  Log("Testing normal_order: spin sorting");
  alg = spin_algebra(2);
  {
    // Sz{1}*Sz{0} -> Sz{0}*Sz{1}
    auto r = normal_order(1.0 * (Op("Sz", 1) * Op("Sz", 0)), alg);
    REQUIRE(isapprox(r, 1.0 * (Op("Sz", 0) * Op("Sz", 1)), alg));
  }
  {
    // S+{1}*S-{0} -> S-{0}*S+{1}  (no sign: spin operators commute)
    auto r = normal_order(1.0 * (Op("S+", 1) * Op("S-", 0)), alg);
    REQUIRE(isapprox(r, 1.0 * (Op("S-", 0) * Op("S+", 1)), alg));
  }

  // -------------------------------------------------------------------------
  // ScalarChirality expansion
  // -------------------------------------------------------------------------
  Log("Testing normal_order: ScalarChirality");
  alg = spin_algebra(3);
  {
    // (i/2)*[S+{0}*S-{1}*Sz{2} - S-{0}*S+{1}*Sz{2} + ...]
    auto r = normal_order(OpSum(Op("ScalarChirality", {0, 1, 2})), alg);
    complex I(0.0, 1.0);
    auto expected = (I * 0.5) * (Op("S+", 0) * Op("S-", 1) * Op("Sz", 2)) +
                    (-I * 0.5) * (Op("S-", 0) * Op("S+", 1) * Op("Sz", 2)) +
                    (I * 0.5) * (Op("Sz", 0) * Op("S+", 1) * Op("S-", 2)) +
                    (-I * 0.5) * (Op("Sz", 0) * Op("S-", 1) * Op("S+", 2)) +
                    (I * 0.5) * (Op("S-", 0) * Op("Sz", 1) * Op("S+", 2)) +
                    (-I * 0.5) * (Op("S+", 0) * Op("Sz", 1) * Op("S-", 2));
    REQUIRE(isapprox(r, expected, alg));
  }

  // -------------------------------------------------------------------------
  // Electron algebra — expansion
  // -------------------------------------------------------------------------
  Log("Testing normal_order: electron expansion");
  {
    // Nup{0} -> Cdagup{0}*Cup{0}
    alg = electron_algebra(1);
    auto r = normal_order(OpSum(Op("Nup", 0)), alg);
    REQUIRE(isapprox(r, 1.0 * (Op("Cdagup", 0) * Op("Cup", 0)), alg));
  }
  {
    // Hopup{0,1}: after sorting -Cdagup{1}*Cup{0} -> +Cup{0}*Cdagup{1}
    alg = electron_algebra(2);
    auto r = normal_order(OpSum(Op("Hopup", {0, 1})), alg);
    auto expected = -1.0 * (Op("Cdagup", 0) * Op("Cup", 1)) +
                    1.0 * (Op("Cup", 0) * Op("Cdagup", 1));
    REQUIRE(isapprox(r, expected, alg));
  }

  // -------------------------------------------------------------------------
  // Electron algebra — same-site CAR
  // -------------------------------------------------------------------------
  Log("Testing normal_order: electron CAR");
  alg = electron_algebra(1);
  {
    // Cdagup{0}*Cdagup{0} -> 0
    auto r = normal_order(1.0 * (Op("Cdagup", 0) * Op("Cdagup", 0)), alg);
    REQUIRE(isapprox(r, OpSum{}, alg));
  }
  {
    // Cup{0}*Cdagup{0} -> Id - Cdagup{0}*Cup{0}
    auto r = normal_order(1.0 * (Op("Cup", 0) * Op("Cdagup", 0)), alg);
    REQUIRE(isapprox(r, 1.0 * Op("Id") - 1.0 * (Op("Cdagup", 0) * Op("Cup", 0)),
                     alg));
  }
  {
    // Cdagup{0}*Cdagdn{0} -> -Cdagdn{0}*Cdagup{0}
    auto r = normal_order(1.0 * (Op("Cdagup", 0) * Op("Cdagdn", 0)), alg);
    REQUIRE(isapprox(r, -1.0 * (Op("Cdagdn", 0) * Op("Cdagup", 0)), alg));
  }

  // -------------------------------------------------------------------------
  // Electron algebra — fermionic sorting (-1 per transposition)
  // -------------------------------------------------------------------------
  Log("Testing normal_order: fermionic sorting");
  alg = electron_algebra(2);
  {
    // Cdagup{1}*Cdagup{0} -> -Cdagup{0}*Cdagup{1}
    auto r = normal_order(1.0 * (Op("Cdagup", 1) * Op("Cdagup", 0)), alg);
    REQUIRE(isapprox(r, -1.0 * (Op("Cdagup", 0) * Op("Cdagup", 1)), alg));
  }
  {
    // Cdagup{1}*Cup{0} -> -Cup{0}*Cdagup{1}
    auto r = normal_order(1.0 * (Op("Cdagup", 1) * Op("Cup", 0)), alg);
    REQUIRE(isapprox(r, -1.0 * (Op("Cup", 0) * Op("Cdagup", 1)), alg));
  }

  // -------------------------------------------------------------------------
  // tJ algebra — projected CAR at same site
  // -------------------------------------------------------------------------
  Log("Testing normal_order: tJ projected CAR");
  alg = tj_algebra(1);
  {
    // Cup{0}*Cdagup{0} in tJ -> Id - Cdagup{0}*Cup{0} - Cdagdn{0}*Cdn{0}
    auto r = normal_order(1.0 * (Op("Cup", 0) * Op("Cdagup", 0)), alg);
    auto expected = 1.0 * Op("Id") - 1.0 * (Op("Cdagup", 0) * Op("Cup", 0)) -
                    1.0 * (Op("Cdagdn", 0) * Op("Cdn", 0));
    REQUIRE(isapprox(r, expected, alg));
  }
  {
    // Cdn{0}*Cdagup{0} in tJ -> 0  (no double occupancy)
    auto r = normal_order(1.0 * (Op("Cdn", 0) * Op("Cdagup", 0)), alg);
    REQUIRE(isapprox(r, OpSum{}, alg));
  }
  {
    // Cup{0}*Cdagdn{0} in tJ -> 0  (no double occupancy)
    auto r = normal_order(1.0 * (Op("Cup", 0) * Op("Cdagdn", 0)), alg);
    REQUIRE(isapprox(r, OpSum{}, alg));
  }

  // -------------------------------------------------------------------------
  // HubbardU expansion
  // -------------------------------------------------------------------------
  Log("Testing normal_order: HubbardU / Nupdn");
  alg = electron_algebra(2);
  {
    // The site-free "HubbardU" represents sum_i Nupdn{i}; the electron algebra
    // knows nsites, so it expands it (and then down to elementary operators).
    // It must agree with the explicit sum over sites and leave no HubbardU type.
    auto r = normal_order(OpSum(Op("HubbardU")), alg);
    auto expected =
        normal_order(OpSum(Op("Nupdn", 0)) + OpSum(Op("Nupdn", 1)), alg);
    REQUIRE(isapprox(r, expected, alg));
    for (auto const &[c, mono] : r) {
      for (auto const &op : mono) {
        REQUIRE(op.type() != "HubbardU");
      }
    }
  }
  {
    // Nupdn{0} -> Cdagup{0}*Cup{0}*Cdagdn{0}*Cdn{0}: a single 4-operator
    // monomial on site 0, with no compound operator remaining.
    alg = electron_algebra(1);
    auto r = normal_order(OpSum(Op("Nupdn", 0)), alg);
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
    alg = spin_algebra(2);
    REQUIRE_THROWS(normal_order(OpSum(Op("Sz", 5)), alg));
  }

  Log("Done testing normal_order");
} catch (xdiag::Error const &e) {
  error_trace(e);
}
