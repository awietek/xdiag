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

#ifndef HYDRA_SYMMETRIES_CHARACTERTABLE_
#define HYDRA_SYMMETRIES_CHARACTERTABLE_

#include <string>
#include <complex>
#include <vector>
#include <map>

#include <hydra/symmetries/spacegroup.h>

namespace hydra { namespace symmetries {

    class CharacterTable
    {
      using complex = std::complex<double>;
    public:
      CharacterTable(const SpaceGroup& space_group,
		     const std::vector<std::string>& names,
		     const std::vector<std::vector<int>>& allowed_symmetries,
		     const std::vector<std::vector<complex>>& characters);

      std::vector<std::string> names() const;
      SpaceGroup little_group(std::string name) const;
      int n_symmetries(std::string name) const;
      std::vector<int> allowed_symmetries(std::string name) const; 
      complex character(const std::string& name, const int& n_sym) const;
      std::vector<complex> characters(const std::string& name) const;
      bool is_real(const std::string& name) const;

    private:
      SpaceGroup space_group_;
      std::vector<std::string> names_;
      int n_symmetries_total_;
      std::map<std::string, int> n_symmetries_;
      std::map<std::string, std::vector<int>> allowed_symmetries_;
      std::map<std::string, std::vector<complex>> characters_;
      std::map<std::string, bool> is_real_;
    };

    void Print(const CharacterTable& table);
    CharacterTable read_charactertable(std::string filename);
  }
}

#endif
