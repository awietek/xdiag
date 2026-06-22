// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <vector>

#include <tests/catch.hpp>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace xdiag::algebra;

// Classify an elementary Cdag/C operator into the creation-major sort key
// (creation: Cdag=0 < C=1 ; sector: up=0 < dn=1 ; site). Returns {2,0,0} for a
// non-Cdag/C operator so callers can ignore it.
static std::array<int64_t, 3> no_key(Op const &op) {
  std::string t = op.type();
  if (t == "Cdag" || t == "Cdagup" || t == "Cdagdn") {
    int64_t sec = (t == "Cdagdn") ? 1 : 0;
    return {0, sec, op[0]};
  }
  if (t == "C" || t == "Cup" || t == "Cdn") {
    int64_t sec = (t == "Cdn") ? 1 : 0;
    return {1, sec, op[0]};
  }
  return {2, 0, 0};
}

// Assert every monomial of a normal-ordered OpSum is in creation-major order:
// all creation operators left of all annihilation operators, up before dn,
// ascending site. (Strictly increasing keys ⇒ also no same-mode duplicates.)
static void require_creation_major(OpSum const &compiled) {
  for (auto const &[c, mono] : compiled) {
    (void)c;
    std::array<int64_t, 3> prev = {-1, -1, -1};
    bool started = false;
    for (int64_t k = 0; k < mono.size(); ++k) {
      std::array<int64_t, 3> key = no_key(mono[k]);
      if (key[0] == 2) {
        continue; // not an elementary Cdag/C operator (e.g. Id)
      }
      if (started) {
        REQUIRE(prev < key);
      }
      prev = key;
      started = true;
    }
  }
}

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

// The implementation algebras (the ones feeding the matrix kernels) must bring
// mixed Cdag/C strings into creation-major normal order.
TEST_CASE("normal_order_creation_major", "[operators]") try {
  Log("Testing normal_order: creation-major (implementation algebras)");

  // Electron: a deliberately out-of-order, mixed-sector 4-operator string.
  {
    Algebra alg = electron_implementation_algebra(4);
    auto r = normal_order(1.0 * (Op("Cdn", 1) * Op("Cdagup", 2) * Op("Cup", 0) *
                                 Op("Cdagdn", 3)),
                          alg);
    require_creation_major(r);
  }
  // tJ: an out-of-order string (different sites, both sectors).
  {
    Algebra alg = tj_implementation_algebra(4);
    auto r = normal_order(1.0 * (Op("Cdn", 0) * Op("Cdagup", 1) * Op("Cup", 2) *
                                 Op("Cdagdn", 3)),
                          alg);
    require_creation_major(r);
  }
  // tJ: the projected S- pair stays Cdagdn*Cup (already creation-major).
  {
    Algebra alg = tj_implementation_algebra(2);
    auto r = normal_order(1.0 * (Op("Cdagdn", 1) * Op("Cup", 0)), alg);
    require_creation_major(r);
  }
  // Fermion (single sector): all Cdag before all C.
  {
    Algebra alg = fermion_implementation_algebra(4);
    auto r =
        normal_order(1.0 * (Op("C", 1) * Op("Cdag", 3) * Op("Cdag", 0)), alg);
    require_creation_major(r);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
}

// The basis must MEAN that normal order: a creation-major monomial applied to
// the vacuum equals (a) the product state with the corresponding occupation and
// (b) the same operators applied one at a time (right to left). Checked on the
// number-non-conserving block so every step stays in one block.
template <typename block_t>
static void test_basis_meaning(block_t const &block,
                               std::vector<int64_t> const &target_local_states,
                               Monomial const &mono) {
  int64_t n = block.nsites();
  State vac = product_state(block, std::vector<int64_t>(n, 0));
  State target = product_state(block, target_local_states);

  // (b) the whole creation-major monomial applied to the vacuum
  State whole = apply(mono, vac);

  // (c) the single operators applied iteratively, right to left
  State iter = vac;
  for (int64_t k = mono.size() - 1; k >= 0; --k) {
    iter = apply(mono[k], iter);
  }

  REQUIRE(isapprox(whole, iter));   // composition: whole == iterative
  REQUIRE(isapprox(whole, target)); // and both build the product state
}

TEST_CASE("normal_order_basis_meaning", "[operators]") try {
  Log("Testing normal_order: basis meaning (3-way state agreement)");

  // local states: empty=0, up=1, dn=2 (electron/tJ); occupied=1 (fermion).
  // Target occupation up@0, up@2, dn@1; creation-major monomial
  // Cdagup_0 Cdagup_2 Cdagdn_1.
  Monomial mono_eltj{Op("Cdagup", 0), Op("Cdagup", 2), Op("Cdagdn", 1)};
  test_basis_meaning(Electron(4), {1, 2, 1, 0}, mono_eltj);
  test_basis_meaning(tJ(4), {1, 2, 1, 0}, mono_eltj);

  // Fermion: occupation @0,2,3; creation-major monomial Cdag_0 Cdag_2 Cdag_3.
  Monomial mono_f{Op("Cdag", 0), Op("Cdag", 2), Op("Cdag", 3)};
  test_basis_meaning(Fermion(4), {1, 0, 1, 1}, mono_f);
} catch (xdiag::Error const &e) {
  error_trace(e);
}

// TotalN / TotalNup / TotalNdn / TotalSz are site-free convenience operators
// that expand to the sum of their local operator over all sites.
TEST_CASE("total_operators", "[operators]") try {
  Log("Testing normal_order: TotalN / TotalNup / TotalNdn / TotalSz");

  auto sum_local = [](int64_t n, std::string const &local) {
    OpSum r;
    for (int64_t i = 0; i < n; ++i) {
      r += Op(local, i);
    }
    return r;
  };

  int64_t n = 3;

  // Algebra-level: the total operator normal-orders to the explicit local sum.
  {
    Algebra alg = spin_algebra(n);
    REQUIRE(isapprox(normal_order(OpSum(Op("TotalSz")), alg),
                     normal_order(sum_local(n, "Sz"), alg), alg));
  }
  {
    Algebra alg = electron_algebra(n);
    REQUIRE(isapprox(normal_order(OpSum(Op("TotalNup")), alg),
                     normal_order(sum_local(n, "Nup"), alg), alg));
    REQUIRE(isapprox(normal_order(OpSum(Op("TotalNdn")), alg),
                     normal_order(sum_local(n, "Ndn"), alg), alg));
    REQUIRE(isapprox(normal_order(OpSum(Op("TotalSz")), alg),
                     normal_order(sum_local(n, "Sz"), alg), alg));
  }
  {
    Algebra alg = tj_algebra(n);
    REQUIRE(isapprox(normal_order(OpSum(Op("TotalNup")), alg),
                     normal_order(sum_local(n, "Nup"), alg), alg));
    REQUIRE(isapprox(normal_order(OpSum(Op("TotalNdn")), alg),
                     normal_order(sum_local(n, "Ndn"), alg), alg));
    REQUIRE(isapprox(normal_order(OpSum(Op("TotalSz")), alg),
                     normal_order(sum_local(n, "Sz"), alg), alg));
  }

  // Functional: the assembled matrix equals the explicit local-sum matrix
  // (exercises the implementation algebra + kernels end to end).
  auto check_matrix = [&](auto const &block, std::string const &total,
                          std::string const &local) {
    int64_t ns = block.nsites();
    arma::mat m_total = matrix(OpSum(Op(total)), block);
    arma::mat m_sum = matrix(sum_local(ns, local), block);
    REQUIRE(isapprox(m_total, m_sum));
  };
  check_matrix(Spinhalf(4), "TotalSz", "Sz");
  check_matrix(Electron(4), "TotalSz", "Sz");
  check_matrix(Electron(4), "TotalNup", "Nup");
  check_matrix(Electron(4), "TotalNdn", "Ndn");
  check_matrix(tJ(4), "TotalSz", "Sz");
  check_matrix(tJ(4), "TotalNup", "Nup");
  check_matrix(tJ(4), "TotalNdn", "Ndn");
} catch (xdiag::Error const &e) {
  error_trace(e);
}
