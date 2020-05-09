#include <algorithm>
#include <cassert>

#include "indexsearch.h"

#include <hydra/hilbertspaces/spinhalf.h>

namespace hydra { namespace indexing {

    template <class hilbertspace_t, class index_t>
    IndexSearch<hilbertspace_t, index_t>::IndexSearch
    (const hilbertspace_t& hilbertspace)
      : size_(hilbertspace.size())
    {
      states_.clear();
      state_t previous;
      int ctr = 0;
      for (auto state : hilbertspace)
	{
	  if (ctr++ != 0) assert(state > previous);
	  previous = state;
	  states_.push_back(state);
	}
    }

    template <class hilbertspace_t, class index_t>
    index_t IndexSearch<hilbertspace_t, index_t>::index
    (const state_t& state) const
    { 
      // Do binary search
      return std::lower_bound(states_.begin(), states_.end(), state) - 
	states_.begin(); 
    }

    template class IndexSearch<hilbertspaces::Spinhalf<uint16>, uint16>;
    template class IndexSearch<hilbertspaces::Spinhalf<uint16>, uint32>;
    template class IndexSearch<hilbertspaces::Spinhalf<uint16>, uint64>;

    template class IndexSearch<hilbertspaces::Spinhalf<uint32>, uint16>;
    template class IndexSearch<hilbertspaces::Spinhalf<uint32>, uint32>;
    template class IndexSearch<hilbertspaces::Spinhalf<uint32>, uint64>;

    template class IndexSearch<hilbertspaces::Spinhalf<uint64>, uint16>;
    template class IndexSearch<hilbertspaces::Spinhalf<uint64>, uint32>;
    template class IndexSearch<hilbertspaces::Spinhalf<uint64>, uint64>;

  }  // namespace indexing
}  // namespace hydra
