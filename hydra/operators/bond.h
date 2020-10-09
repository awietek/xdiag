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

#ifndef HYDRA_OPERATORS_BOND_
#define HYDRA_OPERATORS_BOND_

#include <ostream>
#include <string>
#include <vector>

#include <hydra/parameters/parameters.h>

namespace hydra {

/*!
  Datatype for describing a bond and the sites it lives on

  A bond is described by its type (e.g. Heisenberg, Ising, Multispin,
  Hubbard-U et.), it's coupling constant (i.e. the name of the prefactor)
  and the sites it lives on
 */
struct Bond {
  Bond(const std::string &type, const std::string &coupling,
       const std::vector<int> &sites);
  Bond(const std::string &type, const std::string &coupling,
       const std::vector<int> &sites, const parameters::Parameters &parameters);
  inline std::string type() const { return type_; }
  inline std::string coupling() const { return coupling_; }
  inline std::vector<int> sites() const { return sites_; }
  inline int sites(const int &j) const { return sites_[j]; }
  inline int size() const { return (int)sites_.size(); }
  inline bool has_parameters() const { return has_parameters_; }
  inline parameters::Parameters parameters() const { return parameters_; }

  std::string type_;
  std::string coupling_;
  std::vector<int> sites_;
  bool has_parameters_;
  parameters::Parameters parameters_;
};

std::vector<int> common_sites(Bond b1, Bond b2);

/// writes bond to stream
std::ostream &operator<<(std::ostream &out, const Bond &bond);

/// Comparision for bonds, equal if type, coupling and sites are equal
bool operator==(const Bond &lhs, const Bond &rhs);

/// Datatype for grouping type and coupling names
struct TypeCoupling {
  TypeCoupling(const std::string &type, const std::string &coupling);
  inline std::string type() const { return type_; };
  inline std::string coupling() const { return coupling_; };

  std::string type_;
  std::string coupling_;
};

/// writes type_coupling to stream
std::ostream &operator<<(std::ostream &out, const TypeCoupling &tc);

/// Comparision for TypeCoupling, equal if type and coupling are equal
bool operator==(const TypeCoupling &tc1, const TypeCoupling &tc2);

/// Extracts type and coupling from a bond
TypeCoupling type_coupling(const Bond &bond);
} // namespace hydra

#endif
