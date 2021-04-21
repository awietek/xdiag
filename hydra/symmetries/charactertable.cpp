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

#include "charactertable.h"

namespace hydra {

CharacterTable::CharacterTable(
    std::vector<Representation> const &representations)
    : representations_(representations) {

  for (auto rep : representations) {
    auto name = rep.name;
    names_.push_back(name);

    if (rep.symmetries.size() != rep.characters.size())
      HydraLog.err("Error creating CharacterTable "
                   "rep.symmetries.size() != rep.characters.size()");
    n_symmetries_[name] = rep.symmetries.size();
    symmetries_[name] = rep.symmetries;
    characters_[name] = rep.characters;

    // Check if representation is real
    bool real = true;
    std::vector<double> chis_real;
    for (auto chi : rep.characters) {
      if (std::abs(std::imag(chi)) > 1e-8) {
        real = false;
      }
      chis_real.push_back(std::real(chi));
    }
    is_real_[name] = real;
    characters_real_[name] = chis_real;
  }
}

} // namespace hydra
