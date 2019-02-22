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

#include <fstream>
#include <iostream>
#include "couplings.h"

namespace hydra { namespace operators {

    std::vector<std::string> Couplings::couplings() const
    {
      std::vector<std::string> cps;
      for (auto c : couplings_)
	cps.push_back(c.first);
      return cps;
    }


    Couplings read_couplings(std::string filename)
    {
      
      std::map<std::string, complex> coupling_map;

      // Open file and handle error
      std::ifstream File(filename.c_str());
      if(File.fail()) 
	{
	  std::cerr << "Error in read_couplings: " 
		    << "Could not open file with filename ["
		    << filename << "] given. Abort." << std::endl;
	  exit(EXIT_FAILURE);
	}

      // Advance to coupling lines
      std::string tobeparsed;
      getline(File, tobeparsed);
      while ((tobeparsed.find("[Couplings]") == std::string::npos) &&
	     (tobeparsed.find("[couplings]") == std::string::npos))
	getline(File, tobeparsed);

      // read lines until '[' is found or else until EOF
      while (std::getline(File, tobeparsed))
	{
	  if ((tobeparsed.find('[') != std::string::npos))
	    break;
	  // Parse line
	  std::string name;
	  complex val;
	  std::stringstream stream(tobeparsed);
	  stream >> name;
	  stream >> val;
	  coupling_map[name] = val;
	  if (!File.good())
	    break;
	} 

      return Couplings(coupling_map);
    }
      
  }  // namespace operators
}  // namespace hydra

