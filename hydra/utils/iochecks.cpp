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

#include "iochecks.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <algorithm>

namespace hydra { namespace utils {
    
    template <class T>
    void check_if_contained_in(const T& elem, const std::vector<T>& vec, 
			       std::string name)
    {
      assert(vec.size() > 0);
      if (std::find(vec.begin(), vec.end(), elem) == vec.end())
	{
	  std::cerr << "Unknown " << name << ": " << elem << " (expected ";
	  for (int i = 0; i < (int)vec.size() - 1; ++i)
	    std::cerr << vec[i] << "/";
	  std::cerr << vec[vec.size() - 1] << ")\n";
	  exit(EXIT_FAILURE);
	}
    }


    void check_if_file_exists (const std::string& filename) {
      std::ifstream f(filename.c_str());
      if (!f.good())
	{
	  std::cerr << "Could not open file with filename ["
		    << filename << "] given. Abort." << std::endl;
	  exit(EXIT_FAILURE);
	}
    }
    
    template
    void check_if_contained_in<std::string>
    (const std::string& elem, const std::vector<std::string>& vec, 
     std::string name);

  }  // namespace utils
}  // namespace hydra

