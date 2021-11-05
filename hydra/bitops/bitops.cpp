#include "bitops.h"

#include <sstream>

namespace hydra::bitops {

// bits_to_string
template <typename bit_t>
std::string bits_to_string(bit_t bits, int n, bool reverse) {
  std::stringstream s;
  for (int i = 0; i < n; ++i)
    s << gbit(bits, i);
  std::string st = s.str();
  return reverse ? std::string(st.rbegin(), st.rend()) : st;
}

template std::string bits_to_string<int16_t>(int16_t, int, bool);
template std::string bits_to_string<int32_t>(int32_t, int, bool);
template std::string bits_to_string<int64_t>(int64_t, int, bool);
  
template std::string bits_to_string<uint16_t>(uint16_t, int, bool);
template std::string bits_to_string<uint32_t>(uint32_t, int, bool);
template std::string bits_to_string<uint64_t>(uint64_t, int, bool);

} // namespace hydra::bitops
