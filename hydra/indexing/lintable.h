#ifndef HYDRA_INDEXING_INDEXTABLE_
#define HYDRA_INDEXING_INDEXTABLE_

#include <vector>

#include <hydra/utils/typedefs.h>

namespace hydra { namespace indexing {

    template <class hilbertspace_type, class index_type=uint32>
    class IndexTable {
    public:
      using hilbertspace_t = hilbertspace_type;
      using index_t = index_type;
      using state_t = typename hilbertspace_t::state_t;
    
      IndexTable() = default;
      IndexTable(const hilbertspace_t& hilbertspace);

      index_t index(const state_t& state) const { return indices_[state]; }
      state_t state(const index_t& index) const { return states_[index]; }  
      index_t size() const { return size_; }

    private:
      index_t size_;
      int right_state_size;
      int left_state_size;
      std::vector<std::vector<index_t>> indices_;
      std::vector<std::vector<state_t>> states_;  
    };

  }  // namespace indexing
}  // namespace hydra

#endif
