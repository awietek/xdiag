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

#include <algorithm>
#include <iostream>
#include <sstream>
#include <cassert>
#include <fstream>

#include "bondlist.h"

namespace hydra { namespace operators {
  
    namespace detail 
    {

      std::vector<TypeCoupling> get_types_couplings
      (const std::vector<Bond>& bonds)
      {
	std::vector<TypeCoupling> types_couplings;

	for (Bond bond : bonds)
	  {
	    TypeCoupling tc = type_coupling(bond);
	    if (std::find(types_couplings.begin(), types_couplings.end(), tc) ==
		types_couplings.end())
	      types_couplings.push_back(tc);
	  }
	return types_couplings;
      }
  
      std::vector<std::string> get_types
      (const std::vector<Bond>& bonds)
      {
	std::vector<std::string> types;
	for (Bond bond : bonds)
	  if (std::find(types.begin(), types.end(), bond.type()) == types.end())
	    types.push_back(bond.type());
	return types;
      }

      std::vector<std::string> get_couplings
      (const std::vector<Bond>& bonds)
      {
	std::vector<std::string> couplings;
	for (Bond bond : bonds)
	  if (std::find(couplings.begin(), couplings.end(), 
			bond.coupling()) == couplings.end())
	    couplings.push_back(bond.coupling());
	return couplings;
      }

      std::vector<Bond> get_bonds_of_type(const std::vector<Bond>& bonds,
					  const std::string& type)
      {
	std::vector<Bond> bonds_return;
	for (Bond bond : bonds)
	  if (bond.type() == type) bonds_return.push_back(bond);
	return bonds_return;
      }
  
      std::vector<Bond> get_bonds_of_coupling(const std::vector<Bond>& bonds,
					      const std::string& coupling)
      {
	std::vector<Bond> bonds_return;
	for (Bond bond : bonds)
	  if (bond.coupling() == coupling) bonds_return.push_back(bond);
	return bonds_return;
      }
  
      std::vector<Bond> get_bonds_of_type_coupling
      (const std::vector<Bond>& bonds, const TypeCoupling& type_coupling)
      {
	std::vector<Bond> bonds_return;
	for (Bond bond : bonds)
	  if ((bond.type() == type_coupling.type()) && 
	      (bond.coupling() == type_coupling.coupling()))
	    bonds_return.push_back(bond);
	return bonds_return;
      }

    }  // namespace detail

    BondList::BondList(const std::vector<Bond>& bonds)
    : bonds_(bonds)
    {}

    int BondList::n_sites() const
    {
      int n_sites = 0;
      for (Bond bond : bonds_)
	for (int site : bond.sites())
	  n_sites = std::max(n_sites, site  + 1);

      return n_sites;
    }
    
    void BondList::operator<<(const Bond& bond)
    { bonds_.push_back(bond); }

    std::vector<std::string> BondList::types() const 
    { return detail::get_types(bonds_); }
    std::vector<std::string> BondList::couplings() const 
    { return detail::get_couplings(bonds_); }
    std::vector<TypeCoupling> BondList::types_couplings() const 
    { return detail::get_types_couplings(bonds_); };

    BondList BondList::bonds_of_type(const std::string& type) const
    { return BondList(detail::get_bonds_of_type(bonds_, type)); }

    BondList BondList::bonds_of_coupling(const std::string& coupling) const
    { return BondList(detail::get_bonds_of_coupling(bonds_, coupling)); }

    BondList BondList::bonds_of_type_coupling
    (const std::string& type, const std::string& coupling) const
    { 
      return BondList(detail::get_bonds_of_type_coupling
		      (bonds_, TypeCoupling(type, coupling))); 
    }


    BondList read_bondlist(std::string filename)
    {
      std::vector<Bond> bonds;

      // Open file and handle error
      std::ifstream File(filename.c_str());
      if(File.fail()) 
	{
	  std::cerr << "Error in read_bondlist: " 
		    << "Could not open file with filename ["
		    << filename << "] given. Abort." << std::endl;
	  exit(EXIT_FAILURE);
	}

      // Advance to interaction lines
      std::string tobeparsed;
      getline(File, tobeparsed);
      while ((tobeparsed.find("[Interactions]") == std::string::npos) &&
	     (tobeparsed.find("[interactions]") == std::string::npos))
	getline(File, tobeparsed);

      // read lines until '[' is found or else until EOF
      while (std::getline(File, tobeparsed))
	{
	  if ((tobeparsed.find('[') != std::string::npos) || !File.good())
	    break;
	  
	  // Parse line
	  std::string type, coupling;
	  std::vector<int> sites;
	  std::stringstream stream(tobeparsed);
	  stream >> type;
	  stream >> coupling;
	  int n;
	  while(stream >> n) sites.push_back(n);
	  bonds.push_back(Bond(type, coupling, sites));
	} 

      return BondList(bonds);
    }
      
  }  // namespace operators
}  // namespace hydra

