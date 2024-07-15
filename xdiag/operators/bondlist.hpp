#pragma once

#include <string>
#include <vector>
#include <xdiag/operators/bond.hpp>

namespace xdiag {

class BondList {
public:
  BondList() = default;
  explicit BondList(std::vector<Bond> const &bonds);

  int64_t size() const;
  bool defined(std::string name) const;
  Coupling &operator[](std::string name);
  Coupling const &operator[](std::string name) const;
  std::vector<std::string> couplings() const;

  bool isreal() const;
  bool isexplicit() const;

  void operator+=(Bond const &bond);
  void operator+=(BondList const &bonds);
  BondList operator+(Bond const &bond) const;
  BondList operator+(BondList const &bonds) const;

  bool operator==(BondList const &other) const;
  bool operator!=(BondList const &other) const;

  // Iterators
  using iterator_t = typename std::vector<Bond>::iterator;
  using const_iterator_t = typename std::vector<Bond>::const_iterator;
  iterator_t begin();
  iterator_t end();
  const_iterator_t begin() const;
  const_iterator_t end() const;
  const_iterator_t cbegin() const;
  const_iterator_t cend() const;

private:
  std::vector<Bond> bonds_;
  std::map<std::string, Coupling> couplings_;
};

BondList bonds_of_type(std::string type, BondList const &bonds);
BondList make_explicit(BondList const &bonds);

// legacy function, can be removed later
BondList read_bondlist(std::string filename);

} // namespace xdiag
