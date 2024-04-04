#pragma once

#include <map>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <xdiag/symmetries/continuous_group.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

class QNum {
public:
  using group_t = std::variant<U1, PermutationGroup>;
  using irrep_t = std::variant<int, Representation>;

  QNum() = default;
  QNum(U1 group, int number);
  QNum(PermutationGroup const &group, Representation const &irrep);

  bool good() const;
  group_t group() const;
  irrep_t irrep() const;

  bool operator==(QNum const &rhs) const;
  bool operator!=(QNum const &rhs) const;
  QNum operator*(QNum const &other) const;

private:
  bool good_;
  group_t group_;
  irrep_t irrep_;
};

class QN {
public:
  using const_iterator = std::map<std::string, QNum>::const_iterator;

  QN() = default;
  QN(std::string label, QNum const &qnum);
  QN(std::vector<std::pair<std::string, QNum>> const &labels_qnums);

  bool good() const;
  bool good(std::string label) const;

  bool operator==(QN const &rhs) const;
  bool operator!=(QN const &rhs) const;
  QN operator*(QN const &other) const;

  const_iterator begin() const;
  const_iterator end() const;

private:
  std::map<std::string, QNum> qnums_;
};

} // namespace xdiag
