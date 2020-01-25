#include "combinatorics.h"
#include <cstdio>

namespace hydra { namespace combinatorics {

    int64 binomial(const int& n, const int& k)
    {
      if( k > n || k < 0) return 0;
      int64 res = 1;
      for (int i = 1; i <= k; i++) res = (res * (n - i + 1)) / i;
      return res;
    }
 
   
    uint64 get_nth_pattern(const uint64& n, const int& n_sites, 
			   const int& n_upspins)
    {
      uint64 state = 0;
      uint64 counter = n;
      for (int n_varying_bits = n_upspins-1; n_varying_bits >= 0; 
	   --n_varying_bits)
	{
	  uint64 n_combinations = 0;
	  for(int n_allowed_pos = n_varying_bits; n_allowed_pos <= n_sites;
	      ++n_allowed_pos)
	    {
	      n_combinations += binomial(n_allowed_pos,n_varying_bits);

	      if( n_combinations > counter)
		{
		  counter -= n_combinations -
		    binomial(n_allowed_pos,n_varying_bits);
		  state |= (uint64(1) << n_allowed_pos);
		  break;
		}
	    }
	}
      return state;
    }

    uint64 get_n_for_pattern(const uint64& pattern, const int& n_sites, 
			     const int& n_upspins)
    {
      uint64 n=0;
      uint64 workpattern = pattern;
      for (int n_varying_bits = n_upspins-1; n_varying_bits >= 0;
	   --n_varying_bits)
	{
	  for (int i=0; i <= n_sites; ++i)
	    {
	      // MSB is at 2^i
	      if ((uint64(1) << (i+1)) > workpattern)
		{
		  n += binomial(i,n_varying_bits+1);
		  workpattern ^= (uint64(1) << i);
		  break;
		}
	    }
	}
      return n;
    }
/*
    uint64 shift_n_bits(uint64 n, int shiftfrom, int shiftby) {
      return ((n << shifby) & ~((1 << shiftfrom) - 1)) + (n | ((1 << shiftfrom)-1))
    }
    state_t consecutive_zeros(state_t v) {
    int c;  // output: c will count v's trailing zero bits,
    if (v)
    {
      v = (v ^ (v - 1)) >> 1;  // Set v's trailing 0s to 1s and zero rest
      for (c = 0; v; c++)
      {
      v >>= 1;
      }
    } else {
      c = CHAR_BIT * sizeof(v);
    }
    }
*/

    template <class state_t> state_t down_hole_to_up(state_t downspins, state_t holes) {
      using utils::gbit;
      state_t downspin_tmp = downspins;
      state_t upspin = 0;
      int bit_i_am_testing = 0;
      while (holes > 0) {
        if (~downspin_tmp & 1) {
          upspin = ((holes & 1) << bit_i_am_testing) | upspin;
          holes = holes >> 1;
      }
      downspin_tmp = downspin_tmp >> 1;
      ++bit_i_am_testing;
      }
      return upspin;
    }
    template <class state_t> state_t up_hole_to_down(state_t upspins, state_t holes) {
      using utils::gbit;
      state_t upspin_tmp=upspins;
      state_t downspin = 0;
      int bit_i_am_testing = 0;
      while (holes > 0) {
        if (~upspin_tmp & 1) {
          downspin = ((holes & 1) << bit_i_am_testing) | downspin;
          holes = holes >> 1;
      }
      upspin_tmp = upspin_tmp = upspin_tmp >> 1;
      ++bit_i_am_testing;
      }
      return downspin;
    }

    template uint32 down_hole_to_up(uint32, uint32);
    template uint64 down_hole_to_up(uint64, uint64);
    template uint32 up_hole_to_down(uint32, uint32);
    template uint64 up_hole_to_down(uint64, uint64);
}

}

