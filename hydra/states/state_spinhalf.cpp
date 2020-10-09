#include "state_spinhalf.h"
#include <sstream>

namespace hydra {

template <typename bit_t>
std::string String(state_spinhalf<bit_t> state, number_t n_sites) {
  std::stringstream s;
  for (int i = 0; i < n_sites; ++i)
    s << utils::gbit(state.spins, i);
  std::string st = s.str();
  return std::string(st.rbegin(), st.rend());
}

template std::string String<uint16>(state_spinhalf<uint16> state, int n_sites);
template std::string String<uint32>(state_spinhalf<uint32> state, int n_sites);
template std::string String<uint64>(state_spinhalf<uint64> state, int n_sites);

} // namespace hydra
