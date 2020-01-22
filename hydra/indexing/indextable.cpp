#include "indextable.h"
#include "spinhalf.h"

namespace hydra { namespace indexing {

    template <class hilbertspace_t, class index_t>
    IndexTable<hilbertspace_t, index_t>::IndexTable
    (const hilbertspace_t& hilbertspace)
      : size_(hilbertspace.size()),
	indices_(hilbertspace.rawsize(), -1)
    {
      states_.clear();
      index_t idx = 0;

      for (auto state : hilbertspace)
	{
	  states_.push_back(state);
	  indices_[state] = idx;
	  ++idx;
	}
    }

    template class IndexTable<hilbertspaces::Spinhalf<uint16>, uint16>;
    template class IndexTable<hilbertspaces::Spinhalf<uint16>, uint32>;
    template class IndexTable<hilbertspaces::Spinhalf<uint16>, uint64>;

    template class IndexTable<hilbertspaces::Spinhalf<uint32>, uint16>;
    template class IndexTable<hilbertspaces::Spinhalf<uint32>, uint32>;
    template class IndexTable<hilbertspaces::Spinhalf<uint32>, uint64>;

    template class IndexTable<hilbertspaces::Spinhalf<uint64>, uint16>;
    template class IndexTable<hilbertspaces::Spinhalf<uint64>, uint32>;
    template class IndexTable<hilbertspaces::Spinhalf<uint64>, uint64>;

  }  // namespace indexing
}  // namespace hydra
