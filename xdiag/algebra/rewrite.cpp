// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "rewrite.hpp"

#include <xdiag/algebra/collect.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::operators {

// Try all rules against a monomial. Returns the replacement OpSum for
// the monomial (without the outer coeff), or nullopt if nothing matched.
static std::optional<OpSum>
apply_rules_to_mono(Monomial const &mono,
                    std::vector<MonomialRule> const &mrules,
                    std::vector<OpRule> const &orules) {
  // 1. Try monomial rules on the full monomial (first match wins)
  for (auto const &mrule : mrules) {
    auto maybe = mrule(mono);
    if (maybe)
      return maybe;
  }

  // 2. Try op rules left-to-right (first matching position wins)
  int64_t n = mono.size();
  for (int64_t i = 0; i < n; ++i) {
    for (auto const &orule : orules) {
      auto maybe = orule(mono[i]);
      if (maybe) {
        // Sandwich replacement between prefix and suffix
        std::vector<Op> pre_ops, suf_ops;
        pre_ops.reserve(i);
        suf_ops.reserve(n - i - 1);
        for (int64_t j = 0; j < i; ++j)
          pre_ops.push_back(mono[j]);
        for (int64_t j = i + 1; j < n; ++j)
          suf_ops.push_back(mono[j]);
        Monomial prefix(pre_ops), suffix(suf_ops);

        OpSum result;
        for (auto const &[c_r, m_r] : *maybe) {
          Monomial new_mono = prefix * m_r * suffix;
          result += OpSum(c_r, new_mono);
        }
        return result;
      }
    }
  }
  return std::nullopt;
}

OpSum rewrite_once(OpSum const &ops, std::vector<MonomialRule> const &mrules,
                   std::vector<OpRule> const &orules) {
  OpSum result;
  for (auto const &[coeff, mono] : ops) {
    auto maybe = apply_rules_to_mono(mono, mrules, orules);
    if (maybe) {
      // Scale each replacement term by the original coefficient
      for (auto const &[c_r, m_r] : *maybe) {
        result += OpSum(coeff * c_r, m_r);
      }
    } else {
      result += OpSum(coeff, mono);
    }
  }
  return collect(result);
}

OpSum rewrite(OpSum const &ops, std::vector<MonomialRule> const &mrules,
              std::vector<OpRule> const &orules, int64_t max_iter) try {
  OpSum current = ops;
  for (int64_t iter = 0; iter < max_iter; ++iter) {
    OpSum next = rewrite_once(current, mrules, orules);
    if (next == current)
      return next;
    current = std::move(next);
  }
  XDIAG_THROW(fmt::format(
      "rewrite did not reach a fixed point after {} iterations", max_iter));
}
XDIAG_CATCH

// Convenience overloads
OpSum rewrite_once(OpSum const &ops, std::vector<OpRule> const &orules) {
  return rewrite_once(ops, {}, orules);
}
OpSum rewrite_once(OpSum const &ops, std::vector<MonomialRule> const &mrules) {
  return rewrite_once(ops, mrules, {});
}
OpSum rewrite(OpSum const &ops, std::vector<OpRule> const &orules,
              int64_t max_iter) {
  return rewrite(ops, {}, orules, max_iter);
}
OpSum rewrite(OpSum const &ops, std::vector<MonomialRule> const &mrules,
              int64_t max_iter) {
  return rewrite(ops, mrules, {}, max_iter);
}

} // namespace xdiag
