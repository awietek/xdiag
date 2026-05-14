// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/algebra/rewrite.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace xdiag::operators;

// Empty rule vectors with explicit types to avoid overload ambiguity
static std::vector<MonomialRule> no_mrules;
static std::vector<OpRule> no_orules;

TEST_CASE("test_rewrite", "[operators]") try {
  Log("Testing rewrite engine");

  // --- OpRule: simple type rename ---
  {
    // Rule: rename "Sp" -> "Sx" (preserving sites)
    std::vector<OpRule> orules = {[](Op const &op) -> std::optional<OpSum> {
      if (op.type() == "Sp")
        return OpSum(Op("Sx", op.sites()));
      return std::nullopt;
    }};

    OpSum ops = 2.0 * Op("Sp", 0) + 3.0 * Op("Sp", 1);
    OpSum expected = 2.0 * Op("Sx", 0) + 3.0 * Op("Sx", 1);
    REQUIRE(rewrite(ops, no_mrules, orules) == expected);
  }

  // --- OpRule: multi-term expansion, outer coefficient distributes ---
  {
    // Rule: Sp(i) -> Sx(i) + i*Sy(i)
    std::vector<OpRule> orules = {[](Op const &op) -> std::optional<OpSum> {
      if (op.type() == "Sp")
        return OpSum(Op("Sx", op.sites())) +
               complex(0.0, 1.0) * Op("Sy", op.sites());
      return std::nullopt;
    }};

    // 3 * Sp(0) -> 3*Sx(0) + 3i*Sy(0)
    OpSum ops = 3.0 * Op("Sp", 0);
    OpSum expected = 3.0 * Op("Sx", 0) + complex(0.0, 3.0) * Op("Sy", 0);
    REQUIRE(rewrite(ops, no_mrules, orules) == expected);
  }

  // --- OpRule inside a multi-op monomial: prefix and suffix preserved ---
  {
    // Rule: X(i) -> 2*X(i)
    std::vector<OpRule> orules = {[](Op const &op) -> std::optional<OpSum> {
      if (op.type() == "X")
        return 2.0 * Op("X", op.sites());
      return std::nullopt;
    }};

    // A(0) * X(1) * B(2)  ->  2 * A(0) * X(1) * B(2)  (one step only)
    Monomial mono{Op("A", 0), Op("X", 1), Op("B", 2)};
    OpSum ops = 1.0 * mono;
    OpSum expected = 2.0 * mono;
    REQUIRE(rewrite_once(ops, no_mrules, orules) == expected);
  }

  // --- OpRule: leftmost match per rewrite_once; fixed point via rewrite ---
  {
    // Rule: X -> Y
    std::vector<OpRule> orules = {[](Op const &op) -> std::optional<OpSum> {
      if (op.type() == "X")
        return OpSum(Op("Y", op.sites()));
      return std::nullopt;
    }};

    // X(0) * X(1): rewrite_once only replaces the leftmost X
    Monomial mono{Op("X", 0), Op("X", 1)};
    OpSum ops = 1.0 * mono;

    REQUIRE(rewrite_once(ops, no_mrules, orules) ==
            1.0 * Monomial{Op("Y", 0), Op("X", 1)});

    // rewrite to fixed point replaces all
    REQUIRE(rewrite(ops, no_mrules, orules) ==
            1.0 * Monomial{Op("Y", 0), Op("Y", 1)});
  }

  // --- MonomialRule: whole-monomial match ---
  {
    // Rule: monomial that is exactly A(0) -> B(0), others untouched
    std::vector<MonomialRule> mrules = {
        [](Monomial const &mono) -> std::optional<OpSum> {
          if (mono.size() == 1 && mono[0] == Op("A", 0))
            return OpSum(Op("B", 0));
          return std::nullopt;
        }};

    // A(0) matches -> B(0); A(1) does not match -> unchanged
    OpSum ops = 2.0 * Op("A", 0) + 1.0 * Op("A", 1);
    OpSum expected = 2.0 * Op("B", 0) + 1.0 * Op("A", 1);
    REQUIRE(rewrite(ops, mrules, no_orules) == expected);
  }

  // --- MonomialRule: idempotency N(0)*N(0) -> N(0), iterated ---
  {
    // Rule: inline sub-sequence search for N(0)*N(0) anywhere in the monomial
    std::vector<MonomialRule> mrules = {
        [](Monomial const &mono) -> std::optional<OpSum> {
          Op N0("N", 0);
          int64_t n = mono.size();
          for (int64_t s = 0; s <= n - 2; ++s) {
            if (mono[s] == N0 && mono[s + 1] == N0) {
              std::vector<Op> pre, suf;
              for (int64_t k = 0; k < s; ++k)
                pre.push_back(mono[k]);
              for (int64_t k = s + 2; k < n; ++k)
                suf.push_back(mono[k]);
              // prefix * N(0) * suffix
              OpSum result;
              for (auto const &[c, m_r] : 1.0 * Op("N", 0))
                result += OpSum(c, Monomial(pre) * m_r * Monomial(suf));
              return result;
            }
          }
          return std::nullopt;
        }};

    // N(0)^4 -> N(0)
    Monomial nnnn{Op("N", 0), Op("N", 0), Op("N", 0), Op("N", 0)};
    REQUIRE(rewrite(1.0 * nnnn, mrules, no_orules) == 1.0 * Op("N", 0));
  }

  // --- MonomialRule: site-aware pattern (Sp(i)*Sm(i) -> Id(i)) ---
  {
    std::vector<MonomialRule> mrules = {
        [](Monomial const &mono) -> std::optional<OpSum> {
          if (mono.size() == 2 && mono[0].type() == "Sp" &&
              mono[1].type() == "Sm" && mono[0].size() == 1 &&
              mono[1].size() == 1 && mono[0][0] == mono[1][0])
            return OpSum(Op("Id", mono[0][0]));
          return std::nullopt;
        }};

    // Sp(0)*Sm(0) -> Id(0)
    REQUIRE(rewrite(1.0 * Monomial{Op("Sp", 0), Op("Sm", 0)}, mrules,
                    no_orules) == 1.0 * Op("Id", 0));

    // Sp(0)*Sm(1): different sites, no match — passes through unchanged
    OpSum ops2 = 1.0 * Monomial{Op("Sp", 0), Op("Sm", 1)};
    REQUIRE(rewrite(ops2, mrules, no_orules) == ops2);
  }

  // --- MonomialRule takes priority over OpRule ---
  {
    // MonomialRule: X(0) alone -> Z(0)
    std::vector<MonomialRule> mrules = {
        [](Monomial const &mono) -> std::optional<OpSum> {
          if (mono.size() == 1 && mono[0] == Op("X", 0))
            return OpSum(Op("Z", 0));
          return std::nullopt;
        }};

    // OpRule: X(i) -> Y(i)  (would fire if MonomialRule didn't match)
    std::vector<OpRule> orules = {[](Op const &op) -> std::optional<OpSum> {
      if (op.type() == "X")
        return OpSum(Op("Y", op.sites()));
      return std::nullopt;
    }};

    // MonomialRule fires first -> Z(0), not Y(0)
    REQUIRE(rewrite(1.0 * Op("X", 0), mrules, orules) == 1.0 * Op("Z", 0));
  }

  // --- No matching rules: OpSum passes through unchanged ---
  {
    OpSum ops = 2.0 * Op("A", 0) + 3.0 * Op("B", 1);
    REQUIRE(rewrite_once(ops, no_mrules, no_orules) == ops);
    REQUIRE(rewrite(ops, no_mrules, no_orules) == ops);
  }

  // --- Convenience overloads (OpRule-only and MonomialRule-only) ---
  {
    std::vector<OpRule> orules = {[](Op const &op) -> std::optional<OpSum> {
      if (op.type() == "A")
        return OpSum(Op("B", op.sites()));
      return std::nullopt;
    }};

    OpSum ops = 1.0 * Op("A", 0);
    REQUIRE(rewrite_once(ops, orules) == 1.0 * Op("B", 0));
    REQUIRE(rewrite(ops, orules) == 1.0 * Op("B", 0));

    std::vector<MonomialRule> mrules = {
        [](Monomial const &mono) -> std::optional<OpSum> {
          if (mono.size() == 1 && mono[0].type() == "X")
            return OpSum(Op("Y", mono[0].sites()));
          return std::nullopt;
        }};

    OpSum ops2 = 1.0 * Op("X", 0);
    REQUIRE(rewrite_once(ops2, mrules) == 1.0 * Op("Y", 0));
    REQUIRE(rewrite(ops2, mrules) == 1.0 * Op("Y", 0));
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
}
