#ifndef HYDRA_INDEXING_INDEXSEARCH_
#define HYDRA_INDEXING_INDEXSEARCH_

#include <vector>

#include <hydra/utils/typedefs.h>

namespace hydra { namespace indexing {

    template <class hilbertspace_type, class index_type=uint32>
    class IndexSearch {
    public:
      using hilbertspace_t = hilbertspace_type;
      using index_t = index_type;
      using state_t = typename hilbertspace_t::state_t;
  
      IndexSearch() = default;
      IndexSearch(const hilbertspace_t& hilbertspace);

      index_t index(const state_t& state) const;
      state_t state(const index_t& index) const { return states_[index]; }
      index_t size() const { return size_; }
  
    private:
      index_t size_;
      std::vector<state_t> states_;  
    };

  }  // namespace indexing
}  // namespace hydra

#endif
