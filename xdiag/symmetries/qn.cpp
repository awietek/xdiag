#include "qn.hpp"

#include <iterator>
#include <set>

#include <xdiag/utils/error.hpp>

namespace xdiag {

QNum::QNum(U1 group, int number) : good_(true), group_(group), irrep_(number) {}
QNum::QNum(PermutationGroup const &group, Representation const &irrep)
    : good_(true), group_(group), irrep_(irrep) {}

bool QNum::good() const { return good_; }
QNum::group_t QNum::group() const { return group_; }
QNum::irrep_t QNum::irrep() const { return irrep_; }

bool QNum::operator==(QNum const &rhs) const {
  return (rhs.group_ == group_) && (rhs.irrep_ == irrep_);
}
bool QNum::operator!=(QNum const &rhs) const { return !operator==(rhs); }

QNum QNum::operator*(QNum const &other) const {
  return std::visit(
      overload{
          [&](U1 const &g1, U1 const &g2, int const &i1, int const &i2) {
            if (g1 != g2) {
              throw std::runtime_error(
                  "Error multiplying QNum: U(1) groups not identical");
            }
            return QNum(g1, i1 + i2);
          },
          [&](PermutationGroup const &g1, PermutationGroup const &g2,
              Representation const &i1, Representation const &i2) {
            if (g1 != g2) {
              throw std::runtime_error(
                  "Error multiplying QNum: permutation groups not"
                  "identical");
            }

            try {
              return QNum(g1, i1 * i2);
            } catch (...) {
              throw std::runtime_error(
                  "Error multiplying QNum: unable to multiplyRepresentation");
            }
          },
          [&](auto const &g1, auto const &g2, auto const &i1, auto const &i2) {
            (void)g1;
            (void)g2;
            (void)i1;
            (void)i2;
            throw std::runtime_error(
                "Error multiplying QNum: incompatible group or irrep");
            return QNum();
          }},
      group_, other.group_, irrep_, other.irrep_);
}

QN::QN(std::string label, QNum const &qnum) { qnums_[label] = qnum; }
QN::QN(std::vector<std::pair<std::string, QNum>> const &labels_qnums) {
  for (auto [label, qnum] : labels_qnums) {
    qnums_[label] = qnum;
  }
}

bool QN::good() const {
  for (auto [label, qn] : *this) {
    if (!qn.good()) {
      return false;
    }
  }
  return true;
}

bool QN::good(std::string label) const { return qnums_.at(label).good(); }
bool QN::operator==(QN const &rhs) const { return qnums_ == rhs.qnums_; }
bool QN::operator!=(QN const &rhs) const { return !operator==(rhs); }
QN QN::operator*(QN const &other) const {
  std::set<std::string> keys;
  std::transform(begin(), end(), std::inserter(keys, keys.end()),
                 [](auto pair) { return pair.first; });
  std::set<std::string> keys_other;
  std::transform(other.begin(), end(),
                 std::inserter(keys_other, keys_other.end()),
                 [](auto pair) { return pair.first; });
  if (keys != keys_other) {
    throw std::runtime_error(
        "Error multiplying QN: labels are not the same for the two QNs");
  }

  std::vector<std::pair<std::string, QNum>> labels_qnums;
  for (auto [label, qn] : *this) {
    try {
      QNum qnr = qn * other.qnums_.at(label);
      labels_qnums.push_back({label, qnr});
    } catch (...) {
      throw std::runtime_error(
          std::string(
              "Error multiplying QN: unable to multiply QNum with label ") +
          label);
    }
  }
  return QN(labels_qnums);
}
QN::const_iterator QN::begin() const { return qnums_.begin(); }
QN::const_iterator QN::end() const { return qnums_.end(); }

} // namespace xdiag
