#ifndef HYDRA_UTILS_COMBINATORICS_
#define HYDRA_UTILS_COMBINATORICS_

#include "typedefs.h"

namespace hydra { namespace combinatorics {
    
    int64 binomial(const int& n, const int& k);
    
    template <class int_t> 
    inline int_t get_next_pattern(const int_t& v)
    {
      // Bit twiddling Hack from
      // http://graphics.stanford.edu/~seander/bithacks.html
      // #NextBitPermutation
      
      // // Fast version (needs __builtin_ctz(v))
      // int_t t = v | (v - 1); // t gets v's least significant 0
      // return ((t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1)));

      // Slow version (should work everywhere)
      int_t t = (v | (v - 1)) + 1;  
      return v==0 ? ~v : t | ((((t & -t) / (v & -v)) >> 1) - 1);  
    }
    
    uint64 get_nth_pattern(const uint64& n, const int& n_sites, const int& n_upspins);
    uint64 get_n_for_pattern(const uint64& pattern, const int& n_sites, 
			    const int& n_upspins);

  }  // namespace combinatorics
}  // namespace hydra

#endif