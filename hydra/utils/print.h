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

#ifndef HYDRA_UTILS_PRINT_H_
#define HYDRA_UTILS_PRINT_H_

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/spacegroup.h>
#include <hydra/symmetries/charactertable.h>

#define HydraPrint(X) hydra::utils::PrintPretty(#X,X)

namespace hydra { namespace utils {

    inline void PrintPretty(const char* identifier, 
			    const Bond& bond) 
    {
      printf("%s:\n", identifier);
      printf("%s %s ", bond.type().c_str(), bond.coupling().c_str());
      for (auto site : bond.sites())
	printf("%d ", site);
      printf("\n");      
    }

    inline void PrintPretty(const char* identifier, 
			    const BondList& bondlist) 
    {
      printf("%s:\n", identifier);
      for (auto bond : bondlist)
	{
	  printf("%s %s ", bond.type().c_str(), bond.coupling().c_str());
	  for (auto site : bond.sites())
	    printf("%d ", site);
	  printf("\n");  
	}    
    }

    inline void PrintPretty(const char* identifier, 
			    const Couplings& couplings) 
    {
      printf("%s:\n", identifier);
      for (auto coupling : couplings)
	{
	  printf("%s %f %f\n", coupling.first.c_str(), 
		 std::real(coupling.second), std::imag(coupling.second));
	}    
    }

    inline void PrintPretty(const char* identifier, 
			    const SpaceGroup& group)
    {
      int sym_idx=0;
      printf("%s:\n", identifier);
      for (const auto& sym : group.symmetries())
	{
	  printf("[S%d] ", sym_idx);
	  for (auto p : sym) printf("%d ", p);
	  printf("\n");
	  ++sym_idx;
	}
    }

    inline void PrintPretty(const char* identifier, 
			    const CharacterTable& table)
    {
      printf("%s:\n", identifier);
      for (auto name : table.names())
	{
	  printf("[Representation]=%s\n", name.c_str());
	  if ( table.is_real(name) ) printf("REAL\n");
	  else printf("COMPLEX\n");
	  printf("[AllowedOps]=%d\n", table.n_symmetries(name));
	  for (int n_sym : table.allowed_symmetries(name))
	    printf("%d ", n_sym);
	  printf("\n");
	  for (int n_sym=0; n_sym < table.n_symmetries(name); ++n_sym)
	    printf("%f %f\n", std::real(table.character(name, n_sym)),
		   std::imag(table.character(name, n_sym)));
	}
    }

  }
}
#endif
