#include "state_tj.h"
#include "state_spinhalf.h"

namespace hydra {

template <typename bit_t>
std::string String(state_tj<bit_t> state, int n_sites) {
  return String(state_spinhalf<bit_t>({state.ups}), n_sites) + ";" +
         String(state_spinhalf<bit_t>({state.dns}), n_sites);
}

template std::string String<uint16>(state_tj<uint16> state, int n_sites);
template std::string String<uint32>(state_tj<uint32> state, int n_sites);
template std::string String<uint64>(state_tj<uint64> state, int n_sites);

} // namespace hydra
