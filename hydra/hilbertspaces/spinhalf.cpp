#include <math.h>
#include <cassert>
#include <sstream>

#include "spinhalf.h"
#include "bitops.h"

namespace hydra { namespace hilbertspaces {

    template <class state_t>
    Spinhalf<state_t>::Spinhalf(const int& n_sites, const quantumnumber_t& qn)
    : n_sites_(n_sites), n_upspins_(qn)
    { 
      assert(n_upspins_ <= n_sites_);
      assert(n_upspins_ >= 0);
      state_t begin_state = (((state_t)1 << n_upspins_) - 1);
      state_t end_state = begin_state << (n_sites_ - n_upspins_);
      end_state = combinatorics::get_next_pattern<state_t>(end_state);
      begin_ = SpinhalfIterator<state_t>(begin_state);
      end_ = SpinhalfIterator<state_t>(end_state);
    }

    template <class state_t>
    int Spinhalf<state_t>::n_sites() const { return n_sites_; }

    template <class state_t>
    int Spinhalf<state_t>::quantumnumber() const 
    { return n_upspins_; }

    template <class state_t>
    SpinhalfIterator<state_t> Spinhalf<state_t>::begin() const 
    { return begin_; }
    
    template <class state_t>
    SpinhalfIterator<state_t> Spinhalf<state_t>::end() const 
    { return end_; }
    
    template <class state_t>
    uint64 Spinhalf<state_t>::size() const 
    { return combinatorics::binomial(n_sites_, n_upspins_); }

    template <class state_t>
    uint64 Spinhalf<state_t>::rawsize() const 
    { return pow(2, n_sites_); }

    template <class state_t>
    SpinhalfIterator<state_t>::SpinhalfIterator(const state_t& state)
    : current_(state) 
    {}
    
    template <typename state_t>
    std::string PrintSpinhalf(int n_sites, state_t state)
    {
      std::stringstream s;
      for (int i=0; i<n_sites; ++i)
	s << utils::gbit(state, i);
      std::string st = s.str();
      return std::string(st.rbegin(), st.rend());
    }

    template class Spinhalf<uint16>;
    template class Spinhalf<uint32>;
    template class Spinhalf<uint64>;

    template class SpinhalfIterator<uint16>;
    template class SpinhalfIterator<uint32>;
    template class SpinhalfIterator<uint64>;

    template std::string PrintSpinhalf<uint16>(int n_sites, uint16 state);
    template std::string PrintSpinhalf<uint32>(int n_sites, uint32 state);
    template std::string PrintSpinhalf<uint64>(int n_sites, uint64 state);

  }  // namespace hilbertspaces
}  // namespace hydra
