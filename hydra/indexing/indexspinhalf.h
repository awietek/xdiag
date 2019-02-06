#ifndef HYDRA_INDEXING_INDEXSPINHALF_
#define HYDRA_INDEXING_INDEXSPINHALF_

#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/utils/typedefs.h>

namespace hydra { namespace indexing {

    template <class state_type, class index_type=uint32>
    class IndexSpinhalf {
    public:      
      using index_t = index_type;
      using state_t = state_type;
      using hilbertspace_t = hilbertspaces::Spinhalf<state_t>;

    
      IndexSpinhalf() = default;
      IndexSpinhalf(const hilbertspace_t& hilbertspace);

      index_t index(const state_t& state) const;
      state_t state(const index_t& index) const;
      index_t size() const { return size_; }

    private:
      index_t size_;
      int n_sites_;
      int n_upspins_;
    };


  }  // namespace indexing
}  // namespace hydra
#endif
