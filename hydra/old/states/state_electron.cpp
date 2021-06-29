#include "state_electron.h"
#include "state_spinhalf.h"

namespace hydra {

template <typename bit_t>
std::string String(state_electron<bit_t> state, int n_sites) {
  return String(state_spinhalf<bit_t>({state.ups}), n_sites) + ";" +
         String(state_spinhalf<bit_t>({state.dns}), n_sites);
}

template std::string String<uint16>(state_electron<uint16> state, int n_sites);
template std::string String<uint32>(state_electron<uint32> state, int n_sites);
template std::string String<uint64>(state_electron<uint64> state, int n_sites);

} // namespace hydra
