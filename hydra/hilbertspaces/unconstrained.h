#ifndef HYDRA_HILBERTSPACES_UNCONSTRAINED_
#define HYDRA_HILBERTSPACES_UNCONSTRAINED_

#include <string>

#include <hydra/utils/typedefs.h>

namespace hydra { namespace hilbertspaces {
    
    template <int localdim, class state_type=uint64>
    class UnconstrainedIterator {
    public:
      using state_t = state_type;
      
      UnconstrainedIterator() = default;
      UnconstrainedIterator(const state_t& state);      

      inline bool operator==(const UnconstrainedIterator<state_t>& rhs) const
      { return current_ == rhs.current_; } 
      inline bool operator!=(const UnconstrainedIterator<state_t>& rhs) const
      { return !operator==(rhs); }
      inline UnconstrainedIterator& operator++() { 
	current_ = combinatorics::get_next_pattern(current_);
	return *this;
      }
      inline state_t operator*() const { return current_;}

    private:
      state_t current_;
    };

    template <int localdim, class state_type=uint64>
    class Unconstrained {
    public:
      using quantumnumber_t = int;
      using state_t = state_type;
      using iterator_t = UnconstrainedIterator<localdim, state_t>;
  
      Unconstrained(const int& n_sites);

      int n_sites() const;
      quantumnumber_t quantumnumber() const;
      iterator_t begin() const; 
      iterator_t end() const;
      uint64 size() const;
      uint64 rawsize() const;
      
    private:
      int n_sites_;
      quantumnumber_t n_upspins_;
      iterator_t begin_, end_;
    };
    
    template <typename state_t>
    std::string Print(int n_sites, state_t state);

  }  // namespace hilbertspaces
}  // namespace hydra

#endif
