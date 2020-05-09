#include "combinatorics.h"
#include <cstdio>

#include <iostream>
#include <hydra/all.h>

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


    template <class state_t> state_t
    down_hole_to_up(state_t downspins, state_t holes)
    {
      state_t upspins = 0;
      int bit_i_am_testing = 0;
      while (holes)
	{
	  if (~downspins & 1)
	    {
	      upspins |= ((holes & 1) << bit_i_am_testing);
	      holes >>= 1;
	    }
	  downspins >>= 1;
	  ++bit_i_am_testing;
	}
      return upspins;
    }

    template <class state_t>
    state_t up_hole_to_down(state_t upspins, state_t holes)
    {
      state_t downspins = 0;
      int bit_i_am_testing = 0;
      while (holes)
	{
	  if (~upspins & 1)
	    {
	      downspins |= ((holes & 1) << bit_i_am_testing);
	      holes >>= 1;
	    }
	  upspins >>= 1;
	  ++bit_i_am_testing;
	}
      return downspins;
    }

    template <class state_t> 
    state_t up_down_to_hole(state_t ups, state_t downs) 
    {
      state_t holes = 0;
      int bit_i_am_setting = 0;
      while ((ups) || (downs)) 
    	{
    	  if (~ups & 1) 
    	    {
    	      holes |= ((downs & 1) << bit_i_am_setting);
    	      ++bit_i_am_setting;
    	    }
    	  ups >>= 1;
    	  downs >>= 1;
    	}
      return holes;
    }

    template <class state_t> 
    state_t down_up_to_hole(state_t downs, state_t ups) 
    {
      state_t holes = 0;
      int bit_i_am_setting = 0;
      while ((ups) || (downs)) 
    	{
    	  if (~downs & 1) 
    	    {
    	      holes |= ((ups & 1) << bit_i_am_setting);
    	      ++bit_i_am_setting;
    	    }
    	  ups >>= 1;
    	  downs >>= 1;
    	}
      return holes;
    }

    template uint32 down_hole_to_up<uint32>(uint32, uint32);
    template uint64 down_hole_to_up<uint64>(uint64, uint64);
    template uint32 up_hole_to_down<uint32>(uint32, uint32);
    template uint64 up_hole_to_down<uint64>(uint64, uint64);
    template uint32 up_down_to_hole<uint32>(uint32, uint32);
    template uint64 up_down_to_hole<uint64>(uint64, uint64);
    template uint32 down_up_to_hole<uint32>(uint32, uint32);
    template uint64 down_up_to_hole<uint64>(uint64, uint64);

}

}

