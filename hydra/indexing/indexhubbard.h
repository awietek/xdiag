#ifndef HYDRA_INDEXING_INDEXHUBBARD_
#define HYDRA_INDEXING_INDEXHUBBARD_

#include <vector>

#include <hydra/hilbertspaces/hubbard.h>
#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/utils/typedefs.h>

namespace hydra { namespace indexing {

    template <class indexing_t>
    class IndexHubbard {
    public:
      using index_t = typename indexing_t::index_t;
      using single_state_t = typename indexing_t::state_t;
      using state_t = hilbertspaces::hubbard_state<single_state_t>;
      using hilbertspace_t = hilbertspaces::Hubbard<single_state_t>;  
      
      IndexHubbard() = default;
      IndexHubbard(const hilbertspace_t& hilbertspace);

      inline index_t index(const state_t& state) const; 
      inline state_t state(const index_t& index) const;   
      index_t size() const { return size_; }
      
    private:
      index_t size_;
      indexing_t upindexing_, downindexing_;
    };

    template <class indexing_t>
    IndexHubbard<indexing_t>::IndexHubbard
    (const hilbertspace_t& hilbertspace)
      : size_(hilbertspace.size())
    {
      using hilbertspaces::Spinhalf;

      int n_sites = hilbertspace.n_sites();
      Spinhalf<single_state_t> hs_up(n_sites, 
				     hilbertspace.quantumnumber().n_upspins);
      Spinhalf<single_state_t> hs_down(n_sites, 
				       hilbertspace.quantumnumber().n_downspins);
      upindexing_ = indexing_t(hs_up);
      downindexing_ = indexing_t(hs_down);
    }
    
    template <class indexing_t>
    inline typename indexing_t::index_t IndexHubbard<indexing_t>::index
    (const state_t& state) const
    {
      return upindexing_.index(state.upspins)*downindexing_.size() + 
	downindexing_.index(state.downspins);
    }

    template <class indexing_t>
    inline hilbertspaces::hubbard_state<typename indexing_t::state_t> 
    IndexHubbard<indexing_t>::state(const index_t& index) const
    {
      index_t upindex = index / downindexing_.size();
      index_t downindex = index % downindexing_.size();      
      return { upindexing_.state(upindex), downindexing_.state(downindex)};
    }

  }  // namespace indexing
}  // namespace hydra

#endif
