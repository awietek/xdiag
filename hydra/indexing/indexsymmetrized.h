#ifndef HYDRA_INDEXING_INDEXSYMMETRIZED_
#define HYDRA_INDEXING_INDEXSYMMETRIZED_

#include <vector>

#include <hydra/symmetries/charactertable.h>
#include <hydra/utils/typedefs.h>

namespace hydra { namespace indexing {
    using namespace hydra::symmetries;

    template <class hilbertspace_type, class index_type=uint32, bool fermionic=false>
    class IndexSymmetrized {
    public:
      using hilbertspace_t = hilbertspace_type;
      using index_t = index_type;
      using state_t = typename hilbertspace_t::state_t;

    
      IndexSymmetrized() = default;
      IndexSymmetrized(const hilbertspace_t& hilbertspace, 
		       const CharacterTable& character_table,
		       std::string representation_name);

      index_t index(const state_t& state) const;
      std::string representation_name() const { return representation_name_; }
      state_t state(const index_t& index) const { return states_[index]; }  
      double norm(const index_t& index) const { return norms_[index]; } 
      index_t size() const { return (index_t)states_.size(); }
      
      int find_representative(state_t* state) const;
      
    private:
      CharacterTable character_table_;
      std::string representation_name_;
      int n_symmetries_;
      SpaceGroup space_group_;
      std::vector<complex> characters_;

      std::vector<state_t> states_;
      std::vector<double> norms_;  
    };

  }  // namespace indexing
}  // namespace hydra

#endif
