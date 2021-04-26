// Copyright 2018 Alexander Wietek - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef HYDRA_OPERATORS_BONDLIST_
#define HYDRA_OPERATORS_BONDLIST_

#include <ostream>
#include <string>
#include <vector>

#include "bond.h"

namespace hydra {

class BondList {
  using iterator_t = typename std::vector<Bond>::iterator;
  using const_iterator_t = typename std::vector<Bond>::const_iterator;
  friend BondList operator+(BondList const &, BondList const &);

public:
  BondList() = default;
  explicit BondList(std::vector<Bond> const &bonds);
  void operator<<(Bond const &bond);

  int n_sites() const;
  std::vector<std::string> types() const;
  std::vector<std::string> couplings() const;
  std::vector<TypeCoupling> types_couplings() const;
  BondList bonds_of_type(std::string type) const;
  BondList bonds_of_coupling(std::string coupling) const;
  BondList bonds_of_type_coupling(std::string type, std::string coupling) const;

  iterator_t begin() { return bonds_.begin(); }
  iterator_t end() { return bonds_.end(); }
  const_iterator_t begin() const { return bonds_.begin(); }
  const_iterator_t end() const { return bonds_.end(); }
  const_iterator_t cbegin() const { return bonds_.cbegin(); }
  const_iterator_t cend() const { return bonds_.cend(); }

  Bond operator[](int i) const { return bonds_[i]; }
  Bond &operator[](int i) { return bonds_[i]; }
  int size() const { return (int)bonds_.size(); }
  void clear() { bonds_.clear(); }

private:
  std::vector<Bond> bonds_;
};

BondList operator+(BondList const &, BondList const &);
BondList read_bondlist(std::string filename);

} // namespace hydra

#endif
