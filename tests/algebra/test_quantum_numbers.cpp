// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <tests/blocks/electron/testcases_electron.hpp>
#include <tests/catch.hpp>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/algebra/representation.hpp>
#include <xdiag/algebra/symmetrize.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/symmetries/representation_set.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

// raise[action] / lower[action]: the elementary operators carrying charge +1 /
// -1 under a U(1) "action" (nup, ndn, number).
using RaiseLower = std::map<std::string, std::pair<std::string, std::string>>;
// A charge-0 "filler" operator template: {type, number of sites it acts on}.
using Neutral = std::pair<std::string, int>;

// Build a random OpSum in which EVERY monomial carries exactly the charges
// `target`, but with a random number `n_neutral` of charge-0 filler operators
// mixed in (e.g. SzSz{i,j} * S+{k}). |q| raise/lower operators per action set
// the net charge; the fillers and the operator order are randomized. All
// operators sit on distinct sites, so no same-site reductions occur.
static OpSum random_charge_opsum(int64_t nsites,
                                 std::map<std::string, int64_t> const &target,
                                 RaiseLower const &rl,
                                 std::vector<Neutral> const &neutrals,
                                 int n_neutral, int n_mono, std::mt19937 &gen) {
  std::uniform_real_distribution<double> cdist(0.3, 1.0);
  std::uniform_int_distribution<std::size_t> ndist(
      0, neutrals.empty() ? 0 : neutrals.size() - 1);
  OpSum ops;
  for (int m = 0; m < n_mono; ++m) {
    std::vector<int64_t> sites(nsites);
    std::iota(sites.begin(), sites.end(), 0);
    std::shuffle(sites.begin(), sites.end(), gen);
    int64_t idx = 0;
    std::vector<Op> opv;

    // charge-carrying operators
    for (auto const &[action, q] : target) {
      std::string t = (q >= 0) ? rl.at(action).first : rl.at(action).second;
      for (int64_t c = 0; c < std::abs(q); ++c) {
        opv.push_back(Op(t, sites[idx++]));
      }
    }
    // charge-0 fillers (all operators sit on distinct sites, so stop adding
    // fillers once the remaining sites cannot hold the next one)
    for (int e = 0; e < n_neutral && !neutrals.empty(); ++e) {
      Neutral const &nu = neutrals[ndist(gen)];
      if (idx + nu.second > (int64_t)sites.size()) {
        break;
      }
      if (nu.second == 1) {
        opv.push_back(Op(nu.first, sites[idx++]));
      } else {
        opv.push_back(Op(nu.first, {sites[idx], sites[idx + 1]}));
        idx += 2;
      }
    }
    if (opv.empty()) { // charge 0, no fillers: a single charge-0 hop
      auto const &[raise, lower] = rl.begin()->second;
      opv.push_back(Op(raise, sites[idx++]));
      opv.push_back(Op(lower, sites[idx++]));
    }
    std::shuffle(opv.begin(), opv.end(), gen);
    ops += OpSum(Coeff(cdist(gen)), Monomial(opv));
  }
  return ops;
}

static std::map<std::string, int64_t>
random_target(std::vector<std::string> const &actions, int64_t lo, int64_t hi,
              std::mt19937 &gen) {
  std::uniform_int_distribution<int64_t> qd(lo, hi);
  std::map<std::string, int64_t> t;
  for (auto const &a : actions) {
    t[a] = qd(gen);
  }
  return t;
}

// Tests 1 & 2 (U(1) charge): construct random operators -- monomials of varied
// length mixing charge-0 fillers with the charge-carrying operators -- and
// verify representation() recovers the charge, and that the representation of a
// product is the product of the representations.
static void test_charge_rep(std::string const &name, int64_t nsites,
                            std::vector<std::string> const &actions,
                            RaiseLower const &rl,
                            std::vector<Neutral> const &neutrals,
                            algebra::Algebra const &alg, uint32_t seed) {
  Log("quantum numbers: representation (charge) for {}", name);
  std::mt19937 gen(seed);

  for (int n_neutral = 0; n_neutral <= 3; ++n_neutral) {
    for (int trial = 0; trial < 12; ++trial) {
      // Test 1: representation() retrieves the constructed charge.
      auto target = random_target(actions, -1, 1, gen);
      OpSum ops =
          random_charge_opsum(nsites, target, rl, neutrals, n_neutral, 3, gen);
      for (auto const &a : actions) {
        std::optional<Representation> rep =
            algebra::representation(ops, Representation(a, 0), alg);
        REQUIRE(rep);
        REQUIRE(rep->charge() == target.at(a));
      }

      // Test 2: representation(op1 * op2) == representation(op1) *
      // representation(op2).
      auto t1 = random_target(actions, -1, 1, gen);
      auto t2 = random_target(actions, -1, 1, gen);
      OpSum o1 =
          random_charge_opsum(nsites, t1, rl, neutrals, n_neutral, 2, gen);
      OpSum o2 =
          random_charge_opsum(nsites, t2, rl, neutrals, n_neutral, 2, gen);
      OpSum prod = o1 * o2;
      for (auto const &a : actions) {
        auto r1 = algebra::representation(o1, Representation(a, 0), alg);
        auto r2 = algebra::representation(o2, Representation(a, 0), alg);
        auto rp = algebra::representation(prod, Representation(a, 0), alg);
        REQUIRE(r1);
        REQUIRE(r2);
        REQUIRE(rp);
        REQUIRE(isapprox(*rp, (*r1) * (*r2)));
        REQUIRE(rp->charge() == t1.at(a) + t2.at(a));
      }
    }
  }
}

// Test 3 (U(1) charge): apply an operator of definite charge (with fillers) to
// a state in a given sector; the result must live in the sector with the added
// charges.
template <typename MakeBlock>
static void test_charge_apply(std::string const &name, int64_t nsites,
                              std::vector<std::string> const &actions,
                              RaiseLower const &rl,
                              std::vector<Neutral> const &neutrals,
                              std::map<std::string, int64_t> const &q_in,
                              std::map<std::string, int64_t> const &delta,
                              MakeBlock make_block, uint32_t seed) {
  Log("quantum numbers: apply (charge) for {}", name);
  std::mt19937 gen(seed);
  auto block_in = make_block(q_in);
  REQUIRE(dim(block_in) > 0);

  for (int n_neutral = 0; n_neutral <= 3; ++n_neutral) {
    OpSum ops =
        random_charge_opsum(nsites, delta, rl, neutrals, n_neutral, 2, gen);
    auto block_out = block(ops, block_in); // sector arithmetic used by apply
    for (auto const &a : actions) {
      std::optional<int64_t> q = block_out.irreps().charge(a);
      REQUIRE(q);
      REQUIRE(*q == q_in.at(a) + delta.at(a));
    }
    // The applied state must land in exactly that output block.
    State v = random_state(block_in);
    State w = apply(ops, v);
    REQUIRE(isapprox(w.block(), Block(block_out)));
  }
}

TEST_CASE("quantum_numbers_charge", "[algebra]") try {
  // Spinhalf: U(1) "nup" (Sz sector), S+ raises, S- lowers.
  {
    int64_t n = 8;
    RaiseLower rl{{"nup", {"S+", "S-"}}};
    std::vector<Neutral> neut{{"Sz", 1}, {"SzSz", 2}};
    test_charge_rep("Spinhalf", n, {"nup"}, rl, neut,
                    algebra::symmetry_algebra(Spinhalf(n)), 1);
    test_charge_apply("Spinhalf", n, {"nup"}, rl, neut, {{"nup", 4}},
                      {{"nup", 1}},
                      [&](std::map<std::string, int64_t> const &q) {
                        return Spinhalf(n, q.at("nup"));
                      },
                      2);
  }
  // Boson (local dimension 3): U(1) "number", Adag raises, A lowers.
  {
    int64_t n = 6, d = 3;
    RaiseLower rl{{"number", {"Adag", "A"}}};
    std::vector<Neutral> neut{{"N", 1}};
    test_charge_rep("Boson", n, {"number"}, rl, neut,
                    algebra::symmetry_algebra(Boson(n, d)), 3);
    test_charge_apply("Boson", n, {"number"}, rl, neut, {{"number", 3}},
                      {{"number", 1}},
                      [&](std::map<std::string, int64_t> const &q) {
                        return Boson(n, d, q.at("number"));
                      },
                      4);
  }
  // Fermion: U(1) "number", Cdag raises, C lowers.
  {
    int64_t n = 8;
    RaiseLower rl{{"number", {"Cdag", "C"}}};
    std::vector<Neutral> neut{{"N", 1}, {"NN", 2}};
    test_charge_rep("Fermion", n, {"number"}, rl, neut,
                    algebra::symmetry_algebra(Fermion(n)), 5);
    test_charge_apply("Fermion", n, {"number"}, rl, neut, {{"number", 4}},
                      {{"number", -1}},
                      [&](std::map<std::string, int64_t> const &q) {
                        return Fermion(n, q.at("number"));
                      },
                      6);
  }
  // Electron: U(1) "nup" and "ndn".
  {
    int64_t n = 7;
    RaiseLower rl{{"nup", {"Cdagup", "Cup"}}, {"ndn", {"Cdagdn", "Cdn"}}};
    std::vector<Neutral> neut{{"Ntot", 1},  {"Nup", 1},     {"Ndn", 1},
                              {"Sz", 1},     {"SzSz", 2},    {"NtotNtot", 2}};
    test_charge_rep("Electron", n, {"nup", "ndn"}, rl, neut,
                    algebra::symmetry_algebra(Electron(n)), 7);
    test_charge_apply("Electron", n, {"nup", "ndn"}, rl, neut,
                      {{"nup", 3}, {"ndn", 3}}, {{"nup", 1}, {"ndn", -1}},
                      [&](std::map<std::string, int64_t> const &q) {
                        return Electron(n, q.at("nup"), q.at("ndn"));
                      },
                      8);
  }
  // tJ: U(1) "nup" and "ndn" (no double occupancy).
  {
    int64_t n = 8;
    RaiseLower rl{{"nup", {"Cdagup", "Cup"}}, {"ndn", {"Cdagdn", "Cdn"}}};
    std::vector<Neutral> neut{{"Ntot", 1},  {"Nup", 1},     {"Ndn", 1},
                              {"Sz", 1},     {"SzSz", 2},    {"tJSzSz", 2}};
    test_charge_rep("tJ", n, {"nup", "ndn"}, rl, neut,
                    algebra::symmetry_algebra(tJ(n)), 9);
    test_charge_apply("tJ", n, {"nup", "ndn"}, rl, neut,
                      {{"nup", 3}, {"ndn", 3}}, {{"nup", 1}, {"ndn", -1}},
                      [&](std::map<std::string, int64_t> const &q) {
                        return tJ(n, q.at("nup"), q.at("ndn"));
                      },
                      10);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// Combined: an operator carrying BOTH a nontrivial U(1) charge AND a nontrivial
// momentum, applied to a block carrying both. The charge-carrying part is a
// RANDOM charge OpSum (varied-length monomials with charge-0 fillers, all net
// charge `delta`); it is then symmetrize()-d onto a momentum irrep, so the whole
// operator has definite charge `delta` and definite momentum. make_block(q,
// irrep) builds the block in charge sector `q` AND permutation irrep `irrep`.
// This works for every block (elementary and matrix algebra) including
// multi-site monomials and products that wrap around the cyclic group, because
// representation() compares operators with isapprox on the coefficient-weighted
// matrices rather than by exact matrix identity.
template <typename MakeBlock>
static void test_combined(std::string const &name, int64_t nsites,
                          RaiseLower const &rl,
                          std::vector<Neutral> const &neutrals,
                          std::vector<std::string> const &actions,
                          std::map<std::string, int64_t> const &delta,
                          std::map<std::string, int64_t> const &q_in,
                          algebra::Algebra const &alg, MakeBlock make_block,
                          int64_t k_in, uint32_t seed) {
  using xdiag::testcases::electron::get_cyclic_group_irreps_mult;
  Log("quantum numbers: combined charge+permutation for {}", name);
  auto irreps_mults = get_cyclic_group_irreps_mult(nsites);
  std::vector<Representation> const &irreps = std::get<0>(irreps_mults);
  int64_t nk = (int64_t)irreps.size();
  std::mt19937 gen(seed);

  // random charge-`t` OpSum (n_neutral fillers) projected onto momentum irrep k
  auto make_op = [&](std::map<std::string, int64_t> const &t, int64_t k,
                     int n_neutral, int n_mono) {
    return symmetrize(
        random_charge_opsum(nsites, t, rl, neutrals, n_neutral, n_mono, gen),
        irreps[k]);
  };
  // A momentum projection of a random opsum may vanish (on a transitive group
  // all single sites lie in one orbit, so terms can destructively interfere).
  // Such a null projection is not a meaningful test of the quantum number.
  auto nonzero = [&](OpSum const &o) { return normal_order(o, alg).size() > 0; };

  for (int n_neutral = 0; n_neutral <= 1; ++n_neutral) {
    for (int64_t k = 0; k < nk; ++k) {
      OpSum ops = make_op(delta, k, n_neutral, 2);
      if (!nonzero(ops)) {
        continue;
      }

      // Test 1: representation() recovers the charge(s) AND the momentum.
      for (auto const &a : actions) {
        std::optional<Representation> rep =
            algebra::representation(ops, Representation(a, 0), alg);
        REQUIRE(rep);
        REQUIRE(rep->charge() == delta.at(a));
      }
      std::optional<Representation> repk =
          algebra::representation(ops, irreps[k], alg);
      REQUIRE(repk);
      REQUIRE(isapprox(*repk, irreps[k]));

      // Test 3: apply to a block carrying (q_in, momentum k_in); the result must
      // land in the sector (q_in + delta, momentum k_in * k).
      auto block_in = make_block(q_in, irreps[k_in]);
      if (dim(block_in) == 0) {
        continue;
      }
      auto block_out = block(ops, block_in);
      for (auto const &a : actions) {
        std::optional<int64_t> q = block_out.irreps().charge(a);
        REQUIRE(q);
        REQUIRE(*q == q_in.at(a) + delta.at(a));
      }
      std::map<std::string, int64_t> q_out = q_in;
      for (auto const &a : actions) {
        q_out[a] += delta.at(a);
      }
      auto target = make_block(q_out, irreps[k_in] * irreps[k]);
      REQUIRE(isapprox(Block(block_out), Block(target)));

      State v = random_state(block_in);
      State w = apply(ops, v);
      REQUIRE(isapprox(w.block(), Block(block_out)));
    }
  }

  // Test 2: representation(op_k1 * op_k2) multiplies BOTH charge and momentum.
  for (int64_t k1 = 0; k1 < nk; ++k1) {
    for (int64_t k2 = 0; k2 < nk; ++k2) {
      OpSum o1 = make_op(delta, k1, 0, 2);
      OpSum o2 = make_op(delta, k2, 0, 2);
      OpSum prod = o1 * o2;
      if (!nonzero(prod)) {
        continue;
      }
      for (auto const &a : actions) {
        std::optional<Representation> rp =
            algebra::representation(prod, Representation(a, 0), alg);
        REQUIRE(rp);
        REQUIRE(rp->charge() == 2 * delta.at(a));
      }
      std::optional<Representation> rpk =
          algebra::representation(prod, irreps[0], alg);
      REQUIRE(rpk);
      REQUIRE(isapprox(*rpk, irreps[k1] * irreps[k2]));
    }
  }
}

TEST_CASE("quantum_numbers_combined", "[algebra]") try {
  // Spinhalf: S+ carries nup+1; symmetrized random opsum also carries momentum.
  {
    int64_t n = 4;
    RaiseLower rl{{"nup", {"S+", "S-"}}};
    std::vector<Neutral> neut{{"Sz", 1}, {"SzSz", 2}};
    test_combined(
        "Spinhalf", n, rl, neut, {"nup"}, {{"nup", 1}}, {{"nup", 1}},
        algebra::symmetry_algebra(Spinhalf(n)),
        [&](std::map<std::string, int64_t> const &q, Representation const &ir) {
          return Spinhalf(n, q.at("nup"), ir);
        },
        1, 11);
  }
  // Boson: Adag carries number+1 and a momentum.
  {
    int64_t n = 4, d = 3;
    RaiseLower rl{{"number", {"Adag", "A"}}};
    std::vector<Neutral> neut{{"N", 1}};
    test_combined(
        "Boson", n, rl, neut, {"number"}, {{"number", 1}}, {{"number", 1}},
        algebra::symmetry_algebra(Boson(n, d)),
        [&](std::map<std::string, int64_t> const &q, Representation const &ir) {
          return Boson(n, d, q.at("number"), ir);
        },
        1, 12);
  }
  // Electron: Cdagup carries nup+1 (ndn unchanged) and a momentum.
  {
    int64_t n = 4;
    RaiseLower rl{{"nup", {"Cdagup", "Cup"}}, {"ndn", {"Cdagdn", "Cdn"}}};
    std::vector<Neutral> neut{{"Ntot", 1}, {"Sz", 1}};
    test_combined(
        "Electron", n, rl, neut, {"nup", "ndn"}, {{"nup", 1}, {"ndn", 0}},
        {{"nup", 1}, {"ndn", 1}}, algebra::symmetry_algebra(Electron(n)),
        [&](std::map<std::string, int64_t> const &q, Representation const &ir) {
          return Electron(n, q.at("nup"), q.at("ndn"), ir);
        },
        1, 13);
  }
  // tJ: Cdagup carries nup+1 (ndn unchanged) and a momentum.
  {
    int64_t n = 4;
    RaiseLower rl{{"nup", {"Cdagup", "Cup"}}, {"ndn", {"Cdagdn", "Cdn"}}};
    std::vector<Neutral> neut{{"Ntot", 1}, {"Sz", 1}};
    test_combined(
        "tJ", n, rl, neut, {"nup", "ndn"}, {{"nup", 1}, {"ndn", 0}},
        {{"nup", 1}, {"ndn", 1}}, algebra::symmetry_algebra(tJ(n)),
        [&](std::map<std::string, int64_t> const &q, Representation const &ir) {
          return tJ(n, q.at("nup"), q.at("ndn"), ir);
        },
        1, 14);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// Tests 1, 2, 3 (permutation): build operators of definite momentum with
// symmetrize() and check representation() / products / sectors.
template <typename MakeBlock>
static void test_permutation(std::string const &name, int64_t nsites,
                             Op const &op_a, Op const &op_b,
                             algebra::Algebra const &alg, MakeBlock make_block,
                             int64_t k_in) {
  using xdiag::testcases::electron::get_cyclic_group_irreps_mult;
  Log("quantum numbers: permutation for {}", name);
  auto [irreps, mults] = get_cyclic_group_irreps_mult(nsites);
  int64_t nk = (int64_t)irreps.size();

  // Test 1: a symmetrized operator transforms under the irrep it was projected
  // onto.
  for (int64_t k = 0; k < nk; ++k) {
    OpSum ops = symmetrize(op_a, irreps[k]);
    std::optional<Representation> rep =
        algebra::representation(ops, irreps[k], alg);
    REQUIRE(rep);
    REQUIRE(isapprox(*rep, irreps[k]));
  }

  // Test 2: representation(op_k1 * op_k2) == irrep_k1 * irrep_k2.
  for (int64_t k1 = 0; k1 < nk; ++k1) {
    for (int64_t k2 = 0; k2 < nk; ++k2) {
      OpSum o1 = symmetrize(op_a, irreps[k1]);
      OpSum o2 = symmetrize(op_b, irreps[k2]);
      OpSum prod = o1 * o2;
      std::optional<Representation> rp =
          algebra::representation(prod, irreps[0], alg);
      REQUIRE(rp);
      REQUIRE(isapprox(*rp, irreps[k1] * irreps[k2]));
    }
  }

  // Test 3: applying a momentum-k2 operator to a momentum-k1 state yields a
  // momentum-(k1*k2) state.
  for (int64_t k2 = 0; k2 < nk; ++k2) {
    auto block_in = make_block(irreps[k_in]);
    if (dim(block_in) == 0) {
      continue;
    }
    OpSum ops = symmetrize(op_a, irreps[k2]); // op_a is charge-conserving
    auto block_out = block(ops, block_in);
    REQUIRE(isapprox(Block(block_out),
                     Block(make_block(irreps[k_in] * irreps[k2]))));
  }
}

TEST_CASE("quantum_numbers_permutation", "[algebra]") try {
  // Each block uses its production symmetry_algebra(). For Spinhalf / Boson that
  // is the matrix algebra; collect() keeps the scalar coefficient on the term
  // (not folded into the Matrix), so a momentum operator's phase stays on the
  // coefficient and representation() can read it off.

  // Spinhalf chain (Sz conserves nup).
  {
    int64_t n = 4;
    test_permutation("Spinhalf", n, Op("Sz", 0), Op("Sz", 1),
                     algebra::symmetry_algebra(Spinhalf(n)),
                     [&](Representation const &irrep) {
                       return Spinhalf(n, 2, irrep);
                     },
                     1);
  }
  // Boson chain (N conserves number).
  {
    int64_t n = 4, d = 3;
    test_permutation("Boson", n, Op("N", 0), Op("N", 1),
                     algebra::symmetry_algebra(Boson(n, d)),
                     [&](Representation const &irrep) {
                       return Boson(n, d, 2, irrep);
                     },
                     1);
  }
  // Electron chain.
  {
    int64_t n = 4;
    test_permutation("Electron", n, Op("Ntot", 0), Op("Ntot", 1),
                     algebra::symmetry_algebra(Electron(n)),
                     [&](Representation const &irrep) {
                       return Electron(n, 1, 1, irrep);
                     },
                     1);
  }
  // tJ chain.
  {
    int64_t n = 4;
    test_permutation("tJ", n, Op("Ntot", 0), Op("Ntot", 1),
                     algebra::symmetry_algebra(tJ(n)),
                     [&](Representation const &irrep) {
                       return tJ(n, 1, 1, irrep);
                     },
                     1);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
