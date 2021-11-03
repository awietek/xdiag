#include "up_down_hole.h"

namespace hydra::combinatorics {

template <class bit_t> bit_t down_hole_to_up(bit_t downspins, bit_t holes) {
  bit_t upspins = 0;
  int bit_i_am_testing = 0;
  while (holes) {
    if (~downspins & 1) {
      upspins |= ((holes & 1) << bit_i_am_testing);
      holes >>= 1;
    }
    downspins >>= 1;
    ++bit_i_am_testing;
  }
  return upspins;
}

template <class bit_t> bit_t up_hole_to_down(bit_t upspins, bit_t holes) {
  bit_t downspins = 0;
  int bit_i_am_testing = 0;
  while (holes) {
    if (~upspins & 1) {
      downspins |= ((holes & 1) << bit_i_am_testing);
      holes >>= 1;
    }
    upspins >>= 1;
    ++bit_i_am_testing;
  }
  return downspins;
}

template <class bit_t> bit_t up_down_to_hole(bit_t ups, bit_t downs) {
  bit_t holes = 0;
  int bit_i_am_setting = 0;
  while ((ups) || (downs)) {
    if (~ups & 1) {
      holes |= ((downs & 1) << bit_i_am_setting);
      ++bit_i_am_setting;
    }
    ups >>= 1;
    downs >>= 1;
  }
  return holes;
}

template <class bit_t> bit_t down_up_to_hole(bit_t downs, bit_t ups) {
  bit_t holes = 0;
  int bit_i_am_setting = 0;
  while ((ups) || (downs)) {
    if (~downs & 1) {
      holes |= ((ups & 1) << bit_i_am_setting);
      ++bit_i_am_setting;
    }
    ups >>= 1;
    downs >>= 1;
  }
  return holes;
}

template uint16_t down_hole_to_up<uint16_t>(uint16_t, uint16_t);
template uint32_t down_hole_to_up<uint32_t>(uint32_t, uint32_t);
template uint64_t down_hole_to_up<uint64_t>(uint64_t, uint64_t);

template uint16_t up_hole_to_down<uint16_t>(uint16_t, uint16_t);
template uint32_t up_hole_to_down<uint32_t>(uint32_t, uint32_t);
template uint64_t up_hole_to_down<uint64_t>(uint64_t, uint64_t);

template uint16_t up_down_to_hole<uint16_t>(uint16_t, uint16_t);
template uint32_t up_down_to_hole<uint32_t>(uint32_t, uint32_t);
template uint64_t up_down_to_hole<uint64_t>(uint64_t, uint64_t);

template uint16_t down_up_to_hole<uint16_t>(uint16_t, uint16_t);
template uint32_t down_up_to_hole<uint32_t>(uint32_t, uint32_t);
template uint64_t down_up_to_hole<uint64_t>(uint64_t, uint64_t);

} // namespace hydra::combinatorics
