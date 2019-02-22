// Copyright 2019 Alexander Wietek - All Rights Reserved.
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

#ifndef HYDRA_HILBERTSPACES_SPINHALF_
#define HYDRA_HILBERTSPACES_SPINHALF_

#include <string>

#include <hydra/utils/combinatorics.h>
#include <hydra/utils/typedefs.h>

namespace hydra { namespace hilbertspaces {
    
    template <class state_type=uint64>
    class SpinhalfIterator {
    public:
      using state_t = state_type;
      
      SpinhalfIterator() = default;
      SpinhalfIterator(const state_t& state);      

      inline bool operator==(const SpinhalfIterator<state_t>& rhs) const
      { return current_ == rhs.current_; } 
      inline bool operator!=(const SpinhalfIterator<state_t>& rhs) const
      { return !operator==(rhs); }
      inline SpinhalfIterator& operator++() { 
	current_ = combinatorics::get_next_pattern(current_);
	return *this;
      }
      inline state_t operator*() const { return current_;}

    private:
      state_t current_;
    };

    template <class state_type=uint64>
    class Spinhalf {
    public:
      using quantumnumber_t = int;
      using state_t = state_type;
      using iterator_t = SpinhalfIterator<state_t>;
  
      Spinhalf() = default;
      Spinhalf(const int& n_sites, const quantumnumber_t& qn);

      int n_sites() const;
      quantumnumber_t quantumnumber() const;
      iterator_t begin() const; 
      iterator_t end() const;
      uint64 size() const;
      uint64 rawsize() const;
      
    private:
      int n_sites_;
      quantumnumber_t n_upspins_;
      iterator_t begin_, end_;
    };
    
    template <typename state_t>
    std::string PrintSpinhalf(int n_sites, state_t state);

  }  // namespace hilbertspaces
}  // namespace hydra

#endif
