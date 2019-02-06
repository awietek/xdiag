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


  }
}

