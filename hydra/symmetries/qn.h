#pragma once

#include <map>
#include <string>
#include <utility>
#include <vector>
#include <variant>

#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra {

class QNum {
public:
  using group_t = std::variant<std::string, PermutationGroup>;
  using irrep_t = std::variant<int, Representation>;

  QNum() = default;
  QNum(std::string group, int number);
  QNum(PermutationGroup const &group, Representation const &irrep);

  bool operator==(QNum const &rhs) const;
  bool operator!=(QNum const &rhs) const;

private:
  group_t group_;
  irrep_t irrep_;
};

QNum operator*(QNum const &qnum1, QNum const &qnum2);

  
class QN {
public:
  QN() = default;
  QN(std::string label, QNum const& qnum);
  QN(std::vector<std::pair<std::string, QNum>> const& labels_qnums);
  
  bool operator==(QN const &rhs) const;
  bool operator!=(QN const &rhs) const;
private:
  std::vector<std::string> labels_;
  std::map<std::string, QNum> qnums_;
};

} // namespace hydra
