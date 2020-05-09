#include "lintable.h"

#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/utils/combinatorics.h>
#include <hydra/utils/bitops.h>

namespace hydra { namespace indexing {

    using utils::gbits;
    using utils::popcnt;

    template <class hilbertspace_t,  class index_t>
    LinTable<hilbertspace_t, index_t>::LinTable
    (const hilbertspace_t& hilbertspace)
      : total_size_(hilbertspace.size()),

      // This section could be made cleaner if we restrict ourselves to an even
      // number of sites (so the number of left and right sites is always N/2). 
      // It also implicitly assumes that our Hilbert space
      // is spin 1/2, since that's how Lin tables are implemented - in this
      // case, maybe it's best not to use hilbertspace_t as a template variable?
      
      right_size_(pow(2, ceil(hilbertspace.n_sites()/2.0))),
      left_size_(pow(2,hilbertspace.n_sites() - ceil(hilbertspace.n_sites()/2.0))),
      site_divider_(ceil(hilbertspace.n_sites()/2.0)),
      right_indices_(pow(2, ceil(hilbertspace.n_sites()/2.0)), -1), 
      left_indices_(pow(2,hilbertspace.n_sites() - ceil(hilbertspace.n_sites()/2.0)), -1), 
      num_of_right_sites(ceil(hilbertspace.n_sites()/2.0)),
      num_of_left_sites(hilbertspace.n_sites() - ceil(hilbertspace.n_sites()/2.0))
      {
      states_.clear();
      index_t idx = 0;

      // Iterates over left Hilbert space to fill out its indices.
      
      for (state_t left_state=0; left_state < left_size_; left_state++) {

        if (left_state == 0) {
          left_indices_[left_state] = 0;
        } else {
          left_indices_[left_state] = left_indices_[left_state - 1] 
            + hydra::combinatorics::binomial(num_of_right_sites, hilbertspace.quantumnumber() - popcnt(state_t(left_state-1)));
        }

      }
      
      // Iterates through Hilbert space of a given QN to fill out the rest of the indices
      
      for (auto state : hilbertspace)
	    {
        states_.push_back(state);
        right_indices_[gbits(state, site_divider_, 0)] = idx - left_indices_[state >> (site_divider_)];
        ++idx;
      }
      }

    template <class hilbertspace_t,  class index_t>
    index_t LinTable<hilbertspace_t, index_t>::index
    (const state_t& state) const
    {
      return left_indices_[state >> (site_divider_)] +
	right_indices_[gbits(state, site_divider_, 0)];
    }


    template class LinTable<hilbertspaces::Spinhalf<uint16>, uint16>;
    template class LinTable<hilbertspaces::Spinhalf<uint16>, uint32>;
    template class LinTable<hilbertspaces::Spinhalf<uint16>, uint64>;

    template class LinTable<hilbertspaces::Spinhalf<uint32>, uint16>;
    template class LinTable<hilbertspaces::Spinhalf<uint32>, uint32>;
    template class LinTable<hilbertspaces::Spinhalf<uint32>, uint64>;

    template class LinTable<hilbertspaces::Spinhalf<uint64>, uint16>;
    template class LinTable<hilbertspaces::Spinhalf<uint64>, uint32>;
    template class LinTable<hilbertspaces::Spinhalf<uint64>, uint64>;

  }  // namespace indexing
}  // namespace hydra
