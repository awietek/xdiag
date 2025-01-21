#include "isapprox.hpp"

#include <xdiag/operators/logic/order.hpp>
#include <xdiag/utils/matrix.hpp>

namespace xdiag {

bool isapprox(Op const &op1, Op const &op2, double rtol, double atol) try {
  auto [a1, o1] = order(op1);
  auto [a2, o2] = order(op2);

  if (o1.type() == o2.type()) {
    if (o1.hassites() && o2.hassites()) {
      auto s1 = o1.sites();
      auto s2 = o2.sites();
      if (s1 == s2) {
        if (o1.hasmatrix() && o2.hasmatrix()) {
          return isapprox(o1.matrix() * a1, o2.matrix() * a2, rtol, atol);
        } else if (!o1.hasmatrix() && !o2.hasmatrix()) {
          return isapprox(a1, a2, rtol, atol);
        } else {
          return false;
        }
      } else { // sites disagree
        return false;
      }
    } else if (!o1.hassites() && !o2.hassites()) {
      return isapprox(a1, a2, rtol, atol);
    } else { // one has sites, the other doesnt
      return false;
    }
  } else { // type different
    return false;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool isapprox(OpSum const &ops1, OpSum const &ops2, double rtol,
              double atol) try {
  auto t1 = order(ops1).terms();
  auto t2 = order(ops2).terms();
  if (t1.size() != t2.size()) {
    return false;
  } else {
    for (int64_t i = 0; i < (int64_t)t1.size(); ++i) {
      Scalar a1 = t1[i].first.scalar();
      Scalar a2 = t2[i].first.scalar();
      Op o1 = t1[i].second;
      Op o2 = t2[i].second;
      if (!isapprox(a1, a2, rtol, atol) || !isapprox(o1, o2, rtol, atol)) {
        return false;
      }
    }
    return true;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::optional<Scalar> isapprox_multiple(OpSum const &ops1, OpSum const &ops2,
                                        double rtol, double atol) try {
  auto t1 = order(ops1).terms();
  auto t2 = order(ops2).terms();

  if ((t1.size() != t2.size()) || (t1.size() == 0)) {
    return std::nullopt;
  } else {

    // Determine first non-zero coefficient of ops1 and ratio of the terms
    Scalar ratio = 0.;
    int64_t idx0 = 0;
    for (; idx0 < t1.size(); ++idx0) {
      Scalar a01 = t1[idx0].first.scalar();
      if (abs(a01) > 1e-12) {
        Scalar a02 = t2[idx0].first.scalar();
        ratio = a02 / a01;
        break;
      };
    }
    if (idx0 == t1.size()) {
      XDIAG_THROW("All coefficients in first operator are zero.");
    }

    // See whether all other coefficients have the same ratio
    for (int64_t i = 0; i < (int64_t)t1.size(); ++i) {
      Scalar a1 = t1[i].first.scalar();
      Scalar a2 = t2[i].first.scalar();
      Op o1 = t1[i].second;
      Op o2 = t2[i].second;
      if (!isapprox(a1 * ratio, a2, rtol, atol) ||
          !isapprox(o1, o2, rtol, atol)) {
        return std::nullopt;
      }
    }
    return ratio;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
