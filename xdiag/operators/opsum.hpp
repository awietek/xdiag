#pragma once

#include <string>
#include <vector>
#include <xdiag/operators/op.hpp>

namespace xdiag {

class OpSum {
public:
  using value_type = Op;
  
  OpSum() = default;
  explicit OpSum(std::vector<Op> const &ops);

  int64_t size() const;
  bool defined(std::string name) const;
  Coupling &operator[](std::string name);
  Coupling const &operator[](std::string name) const;
  std::vector<std::string> couplings() const;

  bool isreal() const;
  bool isexplicit() const;

  void push_back(Op const &op);
  void operator+=(Op const &op);
  void operator+=(OpSum const &ops);
  OpSum operator+(Op const &op) const;
  OpSum operator+(OpSum const &ops) const;

  bool operator==(OpSum const &other) const;
  bool operator!=(OpSum const &other) const;

  // Iterators

  using iterator_t = typename std::vector<Op>::iterator;
  using const_iterator_t = typename std::vector<Op>::const_iterator;
  iterator_t begin();
  iterator_t end();
  const_iterator_t begin() const;
  const_iterator_t end() const;
  const_iterator_t cbegin() const;
  const_iterator_t cend() const;

private:
  std::vector<Op> ops_;
  std::map<std::string, Coupling> couplings_;
};

OpSum ops_of_type(std::string type, OpSum const &ops);
OpSum make_explicit(OpSum const &ops);

// legacy function, can be removed later
OpSum read_opsum(std::string filename);

} // namespace xdiag
