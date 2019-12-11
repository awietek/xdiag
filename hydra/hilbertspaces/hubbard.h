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

#ifndef HYDRA_HILBERTSPACES_HUBBARD_
#define HYDRA_HILBERTSPACES_HUBBARD_

#include <iostream>
#include <string>
#include <vector>

#include <hydra/hilbertspaces/spinhalf.h>

namespace hydra { namespace hilbertspaces {

    template <class state_type>
    struct hubbard_state {
      using state_t = state_type;
      state_t upspins;
      state_t downspins;
    };

    template <class state_t>
    inline bool operator<(const hubbard_state<state_t>& s1, 
			  const hubbard_state<state_t>& s2)
    { return ((s1.upspins < s2.upspins) ||
	      ((s1.upspins == s2.upspins) && (s1.downspins < s2.downspins))); }

    template <class state_t>
    inline bool operator==(const hubbard_state<state_t>& s1, 
			   const hubbard_state<state_t>& s2)
    { return ((s1.upspins == s2.upspins) && (s1.downspins == s2.downspins)); }

    /*!
      Data structure for quantum numbers of Hubbard-type Hilbert space.
      The quantum number for this type of Hilbert space consists of the number
      of up spins and down spins.
    */
    struct hubbard_qn {
      int n_upspins;
      int n_downspins;
    };

    std::ostream& operator<<(std::ostream& os, const hubbard_qn& qn);

    template <class state_type=uint32>
    class HubbardIterator {
    public:
      using state_t = state_type;
      
      HubbardIterator() = default;
      HubbardIterator(const Spinhalf<state_t> up, 
		      const Spinhalf<state_t> down,
		      const hubbard_state<state_t>& state);      

      inline bool operator==(const HubbardIterator<state_t>& rhs) const
      { return ((up_iter_ == rhs.up_iter_) && (down_iter_ == rhs.down_iter_)); } 
      inline bool operator!=(const HubbardIterator<state_t>& rhs) const
      { return !operator==(rhs); }
      inline HubbardIterator& operator++() { 
	++down_iter_;
	if (down_iter_ == down_end_)
	  {
	    down_iter_ = down_begin_;
	    ++up_iter_;
	  }
	return *this;
      }
      inline hubbard_state<state_t> operator*() const 
      { return {*up_iter_, *down_iter_}; }

    private:
      int n_sites_;
      SpinhalfIterator<state_t> down_begin_, down_end_;
      SpinhalfIterator<state_t> down_iter_, up_iter_;
    };

    /*!
      Class to generate all configurations of a Hubbard-type Hilbert space
      given the number of upspins and downspins on a given number of sites.

      Usage:
      @code
      #include <hydra/hilberspaces/hubbard.h>

      using Hubbard = hydra::hilbertspaces::Hubbard<>;
      using hydra::hilbertspaces::hubbard_qn;
      using hydra::hilbertspaces::Print;

      int n_upspins = 4;
      int n_downspins = 4;
      int n_sites = 8;

      hubbard_qn qn = {n_upspins, n_downspins};
      Hubbard hs(n_sites, qn);
      int ctr=0;
      for (auto state : hs)
        std::cout << Print(n_sites, state) << std::endl;
      
      @endcode
      @tparam state_type elementary datatype to store up and down configurations
                         one of uint16, uint32, uint64
    */
    template <class state_type=uint32>
    class Hubbard {
    public:
      using quantumnumber_t = hubbard_qn;
      using state_t = hubbard_state<state_type>;
      using iterator_t = HubbardIterator<state_type>;
  
      Hubbard(const int& n_sites, const quantumnumber_t& qn);

      /// returns number of sites of the Hilbert space
      int n_sites() const;

      /// returns the quantum number of the Hilbert space
      quantumnumber_t quantumnumber() const;

      /// returns iterator to first configuration of Hilbertspace
      iterator_t begin() const; 
      
      /// returns iterator to next-to-last configuration of Hilbertspace
      iterator_t end() const;

      /// returns dimension of the Hilbert space with given quantum numbers
      uint64 size() const;

      /// returns raw dimension with all quantum numbers, i.e. pow(4, n_sites)
      uint64 rawsize() const;
      
    private:
      int n_sites_;
      quantumnumber_t qn_;
      iterator_t begin_, end_;
    };

    template <typename state_t>
    std::string print(int n_sites, hubbard_state<state_t> state);
      

  }  // namespace hilbertspaces
}  // namespace hydra

#endif
