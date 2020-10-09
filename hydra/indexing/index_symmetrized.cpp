#include <cassert>
#include <limits>
#include <iostream>

#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/hilbertspaces/siteoperations.h>

#include "indexsymmetrized.h"

namespace hydra { namespace indexing {
    
    template <class hilbertspace_t, class index_t, bool fermionic>
    IndexSymmetrized<hilbertspace_t, index_t, fermionic>::IndexSymmetrized
    (const hilbertspace_t& hilbertspace, 
     const CharacterTable& character_table,
     std::string representation_name)
      : character_table_(character_table),
	representation_name_(representation_name),
	n_symmetries_(character_table_.n_symmetries(representation_name)),
	space_group_(character_table_.little_group(representation_name)),
	characters_(character_table_.characters(representation_name))
    {

      for (auto state : hilbertspace)
	{
	  // Find out whether state is a representative
	  state_t representative = state;
	  find_representative(&representative);

	  if (representative == state)
	    {
	      // Compute norm
	      complex amplitude = 0.0;
	      for (int n_sym = 0; n_sym < n_symmetries_; ++n_sym)
		{
		  state_t derivedstate =
		    space_group_.apply<hilbertspace_t>(n_sym, representative);
		  if (derivedstate == representative)
		    {
		      if (fermionic)
			amplitude += 
			  space_group_.fermi_sign<hilbertspace_t>
			  (n_sym, representative) * characters_[n_sym];
		      else
			amplitude += characters_[n_sym];
		    }
		}
	      
	      // Append state if it's a represenative with non-zero norm
	      if (std::abs(amplitude) > 1e-10)
		{
		  states_.push_back(representative);
		  norms_.push_back(std::sqrt(std::abs(amplitude)));
		}

	    }  // if (representative == state)
	}  // loop over states

      assert(norms_.size()==states_.size());

    }

    template <class hilbertspace_t, class index_t, bool fermionic>
    index_t IndexSymmetrized<hilbertspace_t, index_t, fermionic>::index
    (const state_t& state) const 
    { 
      // Do binary search
      auto it = std::lower_bound(states_.begin(), states_.end(), state);
      if (it != states_.end() && !(state < *it))
	return it - states_.begin();
      else
	return -1; 
     }
    
    template <class hilbertspace_t, class index_t, bool fermionic>
    int IndexSymmetrized<hilbertspace_t, index_t,fermionic>::find_representative
    (state_t* state) const
    {
      state_t representative = *state;
      int rep_sym = 0;

      // go through all symmetries to check whether smaller state exists
      for (int n_sym = 0; n_sym < n_symmetries_; ++n_sym)
	{
	  state_t derivedstate =
	    space_group_.apply<hilbertspace_t>(n_sym, *state);
	  if (derivedstate < representative)
	    {
	      representative = derivedstate;
	      rep_sym = n_sym;
	    }
	}    

      assert(representative <= *state);
      *state = representative;
      return rep_sym;
    }


    template 
    class IndexSymmetrized<SpinHalf<uint16>, uint16, false>;
    template 
    class IndexSymmetrized<SpinHalf<uint16>, uint32, false>;
    template 
    class IndexSymmetrized<SpinHalf<uint16>, uint64, false>;
    template 
    class IndexSymmetrized<SpinHalf<uint16>, uint16, true>;
    template 
    class IndexSymmetrized<SpinHalf<uint16>, uint32, true>;
    template 
    class IndexSymmetrized<SpinHalf<uint16>, uint64, true>;

    template 
    class IndexSymmetrized<SpinHalf<uint32>, uint16, false>;
    template 
    class IndexSymmetrized<SpinHalf<uint32>, uint32, false>;
    template 
    class IndexSymmetrized<SpinHalf<uint32>, uint64, false>;
    template 
    class IndexSymmetrized<SpinHalf<uint32>, uint16, true>;
    template 
    class IndexSymmetrized<SpinHalf<uint32>, uint32, true>;
    template 
    class IndexSymmetrized<SpinHalf<uint32>, uint64, true>;

    template 
    class IndexSymmetrized<SpinHalf<uint64>, uint16, false>;
    template 
    class IndexSymmetrized<SpinHalf<uint64>, uint32, false>;
    template 
    class IndexSymmetrized<SpinHalf<uint64>, uint64, false>;
    template 
    class IndexSymmetrized<SpinHalf<uint64>, uint16, true>;
    template 
    class IndexSymmetrized<SpinHalf<uint64>, uint32, true>;
    template 
    class IndexSymmetrized<SpinHalf<uint64>, uint64, true>;

  }  // namespace indexing
}  // namespace hydra
