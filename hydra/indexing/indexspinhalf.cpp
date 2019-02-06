#include "indexspinhalf.h"

#include "combinatorics.h"

namespace hydra { namespace indexing {

    template <class state_t, class index_t>
    IndexSpinhalf<state_t, index_t>::IndexSpinhalf
    (const hilbertspace_t& hilbertspace)
      : size_(hilbertspace.size()),
	n_sites_(hilbertspace.n_sites()),
	n_upspins_(hilbertspace.quantumnumber())
    {}

    template <class state_t, class index_t>
    index_t IndexSpinhalf<state_t, index_t>::index(const state_t& state) const
    { 
      return combinatorics::get_n_for_pattern((uint64)state, n_sites_, 
					      n_upspins_); 
    }

    template <class state_t, class index_t>
    state_t IndexSpinhalf<state_t, index_t>::state(const index_t& index) const
    { 
      return combinatorics::get_nth_pattern((uint64)index, n_sites_, 
					    n_upspins_); 
    }

    template class IndexSpinhalf<uint16, uint16>;
    template class IndexSpinhalf<uint16, uint32>;
    template class IndexSpinhalf<uint16, uint64>;

    template class IndexSpinhalf<uint32, uint16>;
    template class IndexSpinhalf<uint32, uint32>;
    template class IndexSpinhalf<uint32, uint64>;

    template class IndexSpinhalf<uint64, uint16>;
    template class IndexSpinhalf<uint64, uint32>;
    template class IndexSpinhalf<uint64, uint64>;

  }  // namespace indexing
}  // namespace hydra
