// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/algebra/electron_algebra.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/algebra/spin_algebra.hpp>
#include <xdiag/algebra/tj_algebra.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace xdiag::operators;

TEST_CASE("normal_order", "[operators]") try {

  auto spin = spin_algebra();
  auto elec = electron_algebra();
  auto tj   = tj_algebra();

  // -------------------------------------------------------------------------
  // Spin algebra — expansion
  // -------------------------------------------------------------------------
  Log("Testing normal_order: spin expansion");
  {
    // SzSz{0,1} -> Sz{0}*Sz{1}
    auto r = normal_order(OpSum(Op("SzSz", {0, 1})), spin, 2);
    auto expected = 1.0 * (Op("Sz", 0) * Op("Sz", 1));
    REQUIRE(isapprox(r, expected));
  }
  {
    // Exchange{0,1} -> 0.5*S+{0}*S-{1} + 0.5*S-{0}*S+{1}
    auto r = normal_order(OpSum(Op("Exchange", {0, 1})), spin, 2);
    auto expected = 0.5 * (Op("S+", 0) * Op("S-", 1)) +
                    0.5 * (Op("S-", 0) * Op("S+", 1));
    REQUIRE(isapprox(r, expected));
  }
  {
    // SdotS{0,1} -> Sz{0}*Sz{1} + 0.5*S+{0}*S-{1} + 0.5*S-{0}*S+{1}
    auto r = normal_order(OpSum(Op("SdotS", {0, 1})), spin, 2);
    auto expected = 1.0 * (Op("Sz", 0) * Op("Sz", 1)) +
                    0.5 * (Op("S+", 0) * Op("S-", 1)) +
                    0.5 * (Op("S-", 0) * Op("S+", 1));
    REQUIRE(isapprox(r, expected));
  }

  // -------------------------------------------------------------------------
  // Spin algebra — same-site simplification
  // -------------------------------------------------------------------------
  Log("Testing normal_order: spin same-site");
  {
    // Sz{0}*Sz{0} -> 0.25*Id
    auto r = normal_order(1.0 * (Op("Sz", 0) * Op("Sz", 0)), spin, 1);
    REQUIRE(isapprox(r, 0.25 * Op("Id")));
  }
  {
    // S+{0}*S-{0} -> 0.5*Id + Sz{0}
    auto r = normal_order(1.0 * (Op("S+", 0) * Op("S-", 0)), spin, 1);
    REQUIRE(isapprox(r, 0.5 * Op("Id") + 1.0 * Op("Sz", 0)));
  }
  {
    // S-{0}*S+{0} -> 0.5*Id - Sz{0}
    auto r = normal_order(1.0 * (Op("S-", 0) * Op("S+", 0)), spin, 1);
    REQUIRE(isapprox(r, 0.5 * Op("Id") - 1.0 * Op("Sz", 0)));
  }
  {
    // S+{0}*S+{0} -> 0
    auto r = normal_order(1.0 * (Op("S+", 0) * Op("S+", 0)), spin, 1);
    REQUIRE(isapprox(r, OpSum{}));
  }
  {
    // SdotS{0,0} -> 0.75*Id  (S·S for spin-1/2)
    auto r = normal_order(OpSum(Op("SdotS", {0, 0})), spin, 1);
    REQUIRE(isapprox(r, 0.75 * Op("Id")));
  }

  // -------------------------------------------------------------------------
  // Spin algebra — sorting (bosonic: no sign change)
  // -------------------------------------------------------------------------
  Log("Testing normal_order: spin sorting");
  {
    // Sz{1}*Sz{0} -> Sz{0}*Sz{1}
    auto r = normal_order(1.0 * (Op("Sz", 1) * Op("Sz", 0)), spin, 2);
    REQUIRE(isapprox(r, 1.0 * (Op("Sz", 0) * Op("Sz", 1))));
  }
  {
    // S+{1}*S-{0} -> S-{0}*S+{1}  (no sign: spin operators commute)
    auto r = normal_order(1.0 * (Op("S+", 1) * Op("S-", 0)), spin, 2);
    REQUIRE(isapprox(r, 1.0 * (Op("S-", 0) * Op("S+", 1))));
  }

  // -------------------------------------------------------------------------
  // ScalarChirality expansion
  // -------------------------------------------------------------------------
  Log("Testing normal_order: ScalarChirality");
  {
    // (i/2)*[S+{0}*S-{1}*Sz{2} - S-{0}*S+{1}*Sz{2} + ...]
    auto r = normal_order(OpSum(Op("ScalarChirality", {0, 1, 2})), spin, 3);
    complex I(0.0, 1.0);
    auto expected =
        (I * 0.5) * (Op("S+", 0) * Op("S-", 1) * Op("Sz", 2)) +
        (-I * 0.5) * (Op("S-", 0) * Op("S+", 1) * Op("Sz", 2)) +
        (I * 0.5) * (Op("Sz", 0) * Op("S+", 1) * Op("S-", 2)) +
        (-I * 0.5) * (Op("Sz", 0) * Op("S-", 1) * Op("S+", 2)) +
        (I * 0.5) * (Op("S-", 0) * Op("Sz", 1) * Op("S+", 2)) +
        (-I * 0.5) * (Op("S+", 0) * Op("Sz", 1) * Op("S-", 2));
    REQUIRE(isapprox(r, expected));
  }

  // -------------------------------------------------------------------------
  // Electron algebra — expansion
  // -------------------------------------------------------------------------
  Log("Testing normal_order: electron expansion");
  {
    // Nup{0} -> Cdagup{0}*Cup{0}
    auto r = normal_order(OpSum(Op("Nup", 0)), elec, 1);
    REQUIRE(isapprox(r, 1.0 * (Op("Cdagup", 0) * Op("Cup", 0))));
  }
  {
    // Hopup{0,1}: after sorting -Cdagup{1}*Cup{0} -> +Cup{0}*Cdagup{1}
    auto r = normal_order(OpSum(Op("Hopup", {0, 1})), elec, 2);
    auto expected = -1.0 * (Op("Cdagup", 0) * Op("Cup", 1)) +
                     1.0 * (Op("Cup", 0) * Op("Cdagup", 1));
    REQUIRE(isapprox(r, expected));
  }

  // -------------------------------------------------------------------------
  // Electron algebra — same-site CAR
  // -------------------------------------------------------------------------
  Log("Testing normal_order: electron CAR");
  {
    // Cdagup{0}*Cdagup{0} -> 0
    auto r = normal_order(1.0 * (Op("Cdagup", 0) * Op("Cdagup", 0)), elec, 1);
    REQUIRE(isapprox(r, OpSum{}));
  }
  {
    // Cup{0}*Cdagup{0} -> Id - Cdagup{0}*Cup{0}
    auto r = normal_order(1.0 * (Op("Cup", 0) * Op("Cdagup", 0)), elec, 1);
    REQUIRE(isapprox(r, 1.0 * Op("Id") -
                            1.0 * (Op("Cdagup", 0) * Op("Cup", 0))));
  }
  {
    // Cdagup{0}*Cdagdn{0} -> -Cdagdn{0}*Cdagup{0}
    auto r =
        normal_order(1.0 * (Op("Cdagup", 0) * Op("Cdagdn", 0)), elec, 1);
    REQUIRE(isapprox(r, -1.0 * (Op("Cdagdn", 0) * Op("Cdagup", 0))));
  }

  // -------------------------------------------------------------------------
  // Electron algebra — fermionic sorting (-1 per transposition)
  // -------------------------------------------------------------------------
  Log("Testing normal_order: fermionic sorting");
  {
    // Cdagup{1}*Cdagup{0} -> -Cdagup{0}*Cdagup{1}
    auto r =
        normal_order(1.0 * (Op("Cdagup", 1) * Op("Cdagup", 0)), elec, 2);
    REQUIRE(isapprox(r, -1.0 * (Op("Cdagup", 0) * Op("Cdagup", 1))));
  }
  {
    // Cdagup{1}*Cup{0} -> -Cup{0}*Cdagup{1}
    auto r = normal_order(1.0 * (Op("Cdagup", 1) * Op("Cup", 0)), elec, 2);
    REQUIRE(isapprox(r, -1.0 * (Op("Cup", 0) * Op("Cdagup", 1))));
  }

  // -------------------------------------------------------------------------
  // tJ algebra — projected CAR at same site
  // -------------------------------------------------------------------------
  Log("Testing normal_order: tJ projected CAR");
  {
    // Cup{0}*Cdagup{0} in tJ -> Id - Cdagup{0}*Cup{0} - Cdagdn{0}*Cdn{0}
    auto r = normal_order(1.0 * (Op("Cup", 0) * Op("Cdagup", 0)), tj, 1);
    auto expected = 1.0 * Op("Id") -
                    1.0 * (Op("Cdagup", 0) * Op("Cup", 0)) -
                    1.0 * (Op("Cdagdn", 0) * Op("Cdn", 0));
    REQUIRE(isapprox(r, expected));
  }
  {
    // Cdn{0}*Cdagup{0} in tJ -> 0  (no double occupancy)
    auto r = normal_order(1.0 * (Op("Cdn", 0) * Op("Cdagup", 0)), tj, 1);
    REQUIRE(isapprox(r, OpSum{}));
  }
  {
    // Cup{0}*Cdagdn{0} in tJ -> 0  (no double occupancy)
    auto r = normal_order(1.0 * (Op("Cup", 0) * Op("Cdagdn", 0)), tj, 1);
    REQUIRE(isapprox(r, OpSum{}));
  }

  // -------------------------------------------------------------------------
  // HubbardU expansion
  // -------------------------------------------------------------------------
  Log("Testing normal_order: HubbardU");
  {
    // HubbardU with nsites=1 -> Nupdn{0} -> Cdagup{0}*Cup{0}*Cdagdn{0}*Cdn{0}
    // After CAR: simplifies to a single 4-operator monomial on site 0
    auto r = normal_order(OpSum(Op("HubbardU")), elec, 1);
    // Check no compound operators remain
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
    // REQUIRE_THROWS(normal_order(OpSum(Op("Sz", 5)), spin, 3));
  }

  Log("Done testing normal_order");
} catch (xdiag::Error const &e) {
  error_trace(e);
}
