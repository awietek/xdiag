#include <string>
#include <math.h>
#include <cassert>

#include "hubbard.h"

namespace hydra { namespace hilbertspaces {

    template <class state_type>
    Hubbard<state_type>::Hubbard(const int& n_sites, const quantumnumber_t& qn)
    : n_sites_(n_sites), qn_(qn)
    { 
      assert(qn.n_upspins <= n_sites);
      assert(qn.n_upspins >= 0);
      assert(qn.n_downspins <= n_sites);
      assert(qn.n_downspins >= 0);

      Spinhalf<state_type> up(n_sites, qn.n_upspins);
      Spinhalf<state_type> down(n_sites, qn.n_downspins);
      begin_ = HubbardIterator<state_type>(up, down, {*up.begin(), *down.begin()});
      end_ = HubbardIterator<state_type>(up, down, {*up.end(), *down.begin()});
    }
    
    template <class state_t>
    int Hubbard<state_t>::n_sites() const { return n_sites_; }

    std::ostream& operator<<(std::ostream& os, const hubbard_qn& qn)
    { 
      os << "nup: " << qn.n_upspins << ", ndown: " << qn.n_downspins; 
      return os; 
    }

    template <class state_t>
    hubbard_qn Hubbard<state_t>::quantumnumber() const 
    { return qn_; }

    template <class state_t>
    HubbardIterator<state_t> Hubbard<state_t>::begin() const 
    { return begin_; }
    
    template <class state_t>
    HubbardIterator<state_t> Hubbard<state_t>::end() const 
    { return end_; }

    template <class state_t>
    uint64 Hubbard<state_t>::size() const 
    { 
      return combinatorics::binomial(n_sites_, qn_.n_upspins) * 
	combinatorics::binomial(n_sites_, qn_.n_downspins); 
    }

    template <class state_t>
    uint64 Hubbard<state_t>::rawsize() const 
    { return pow(4, n_sites_); }

    template <class state_t>
    HubbardIterator<state_t>::HubbardIterator(const Spinhalf<state_t> up, 
					      const Spinhalf<state_t> down,
					      const hubbard_state<state_t>& state)
      : n_sites_(up.n_sites()),
	down_begin_(down.begin()),
	down_end_(down.end()),
	down_iter_(state.downspins),
	up_iter_(state.upspins)
    { assert(up.n_sites() == down.n_sites()); }




    template class Hubbard<uint16>;
    template class Hubbard<uint32>;
    template class Hubbard<uint64>;

    template class HubbardIterator<uint16>;
    template class HubbardIterator<uint32>;
    template class HubbardIterator<uint64>;

    // template <typename state_t>
    // std::string Print(int n_sites, hubbard_state<state_t> state)
    // { 
    //   return PrintSpinhalf(n_sites, state.upspins) + ";" +
    // 	PrintSpinhalf(n_sites, state.downspins);
    // }
    // template std::string Print<uint16>(int n_sites, 
    // 				       hubbard_state<uint16> state);
    // template std::string Print<uint32>(int n_sites, 
    // 				       hubbard_state<uint32> state);
    // template std::string Print<uint64>(int n_sites, 
    // 				       hubbard_state<uint64> state);


  }  // namespace hilbertspaces
}  // namespace hydra
