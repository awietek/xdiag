#include "lintable_raw.h"

namespace hydra::indexing {

template <class bit_t>
LinTableRaw<bit_t>::LinTableRaw(int n) : n_(n), size_((idx_t)1 << n) {}

template class LinTableRaw<uint16_t>;
template class LinTableRaw<uint32_t>;
template class LinTableRaw<uint64_t>;

} // namespace hydra::indexing
