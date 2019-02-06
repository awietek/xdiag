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

namespace hydra { namespace operators {

    /*!
      BondList is a container for organizing interaction types, coupling 
      coupling constants and the sites, interactions act on.

      Usage:
      \code
      #include <hydra/all.h>
      
      using hydra::operators::BondList;
      using hydra::operators::read_bondlist;

      BondList bondlist = read_bondlist(latticefile);
      int n_sites = bondlist.n_sites();
      BondList hopping_list = bondlist.bonds_of_type("HUBBARDHOP");
      std::vector<std::pair<int, int>> hoppings;
      for (auto bond : hopping_list)
        hoppings.push_back({bond.sites()[0], bond.sites()[1]});
      \endcode
    */
    class BondList
    {
      using iterator_t = typename std::vector<Bond>::iterator;
      using const_iterator_t = typename std::vector<Bond>::const_iterator;

    public:
      /// Create empty Bondlist
      BondList() = default;
      
      /// Create BondList from standard vector of Bond 
      explicit BondList(const std::vector<Bond>& bonds);

      /// Add a Bond to the bondlist
      void operator<<(const Bond& bond);
      
      /// returns on how many sites all the bonds live (i.e. max no. site + 1) 
      int n_sites() const;

      /// returns vector with names of types 
      std::vector<std::string> types() const;

      /// returns vector with names of couplings
      std::vector<std::string> couplings() const;

      /// returns vector names of  with types and couplings
      std::vector<TypeCoupling> types_couplings() const;

      /// returns a BondList with bonds of given type
      BondList bonds_of_type(const std::string& type) const;

      /// returns a BondList with bonds of given coupling
      BondList bonds_of_coupling(const std::string& coupling) const;

      /// returns a BondList with bonds of given type and coupling
      BondList bonds_of_type_coupling
      (const std::string& type, const std::string& coupling) const;

      iterator_t begin() { return bonds_.begin(); }
      iterator_t end() { return bonds_.end(); }
      const_iterator_t begin() const { return bonds_.begin(); }
      const_iterator_t end() const { return bonds_.end(); }
      const_iterator_t cbegin() const { return bonds_.cbegin(); }
      const_iterator_t cend() const { return bonds_.cend(); }
	
    private:
      std::vector<Bond> bonds_;
    };

    /*!
      Reads BondList from filename

      Format should be 
      [Interactions]
      HB J1 0 1
      HB J2 0 2
      ...
      reads a block starting with [Interactions] or [interactions] until next 
      [...] block or if [...] is not found until end of file
    */
    BondList read_bondlist(std::string filename);
      
  }  // namespace operators
}  // namespace hydra

#endif
