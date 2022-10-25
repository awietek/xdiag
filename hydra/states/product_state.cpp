#include "product_state.h"

#include <hydra/indexing/indexing_variants.h>

namespace hydra {

template <typename coeff_t = complex, typename bit_t = std_bit_t>
State<coeff_t, Spinhalf<bit_t>>
product_state(std::vector<std::string> const &local_states,
              Spinhalf<bit_t> const &block) {
  int n_sites = block.n_sites();
  auto state = zero_state(block);

  if (local_states.size() != n_sites) {
    Log.err("Error creating product state: length of local_states does not "
            "match number of sites in block.");
  }

  if (block.symmetric()) {
    Log.err("Error creating product state: cannot create product state on "
            "symmetric block");
  }

  // Create bit pattern
  bit_t bit_pattern = 0;
  for (int s = 0; s < n_sites; ++s) {
    if (local_states[s] == "Up") {
      bit_pattern |= ((bit_t)1 << s);
    } else {
      if (local_states[s] != "Dn") {
        Log.err(
            "Error creating product state: invalid local state encountered: {}",
            local_states[s]);
      }
    }
  }

  // using namespace indexing::spinhalf;
  if (block.sz_conserved()) {
    // auto const &indexing = std::get<IndexingSz<bit_t>>(block.indexing());
    // idx_t idx = indexing.index(bit_pattern);
    // state.vector()(idx) = 1.0;
  } else {
    // auto const &indexing = std::get<IndexingNoSz<bit_t>>(block.indexing());
    // idx_t idx = indexing.index(bit_pattern);
    // state.vector()(idx) = 1.0;
  }

  return state;
}

template State<double, Spinhalf<uint16_t>>
product_state<double, uint16_t>(std::vector<std::string> const &local_states,
                                Spinhalf<uint16_t> const &block);
template State<double, Spinhalf<uint32_t>>
product_state<double, uint32_t>(std::vector<std::string> const &local_states,
                                Spinhalf<uint32_t> const &block);
template State<double, Spinhalf<uint64_t>>
product_state<double, uint64_t>(std::vector<std::string> const &local_states,
                                Spinhalf<uint64_t> const &block);

template State<complex, Spinhalf<uint16_t>>
product_state<complex, uint16_t>(std::vector<std::string> const &local_states,
                                 Spinhalf<uint16_t> const &block);
template State<complex, Spinhalf<uint32_t>>
product_state<complex, uint32_t>(std::vector<std::string> const &local_states,
                                 Spinhalf<uint32_t> const &block);
template State<complex, Spinhalf<uint64_t>>
product_state<complex, uint64_t>(std::vector<std::string> const &local_states,
                                 Spinhalf<uint64_t> const &block);

template <typename coeff_t = complex, typename bit_t = std_bit_t>
State<coeff_t, Electron<bit_t>>
product_state(std::vector<std::string> const &local_states,
              Electron<bit_t> const &block) {
  int n_sites = block.n_sites();
  auto state = zero_state(block);

  if (local_states.size() != n_sites) {
    Log.err("Error creating product state: length of local_states does not "
            "match number of sites in block.");
  }

  if (block.symmetric()) {
    Log.err("Error creating product state: cannot create product state on "
            "symmetric block");
  }

  // Create bit pattern
  bit_t up_pattern = 0;
  bit_t dn_pattern = 0;
  for (int s = 0; s < n_sites; ++s) {
    if (local_states[s] == "Up") {
      up_pattern |= ((bit_t)1 << s);
    } else if (local_states[s] == "Dn") {
      dn_pattern |= ((bit_t)1 << s);
    } else if (local_states[s] == "UpDn") {
      up_pattern |= ((bit_t)1 << s);
      dn_pattern |= ((bit_t)1 << s);
    }
    if (local_states[s] != "Emp") {
      Log.err(
          "Error creating product state: invalid local state encountered: {}",
          local_states[s]);
    }
  }

  // using namespace indexing::electron;
  if (block.charge_conserved()) {
    // auto const &indexing = std::get<IndexingNp<bit_t>>(block.indexing());
    // idx_t idx_up = indexing.index_ups(up_pattern);
    // idx_t idx_dn = indexing.index_dns(dn_pattern);
    // idx_t size_dns = indexing.size_dns();
    // state.vector()(idx_up * size_dns + idx_dn) = 1.0;
  } else {
    // auto const &indexing = std::get<IndexingNoNp<bit_t>>(block.indexing());
    // idx_t idx_up = indexing.index_ups(up_pattern);
    // idx_t idx_dn = indexing.index_dns(dn_pattern);
    // idx_t size_dns = indexing.size_dns();
    // state.vector()(idx_up * size_dns + idx_dn) = 1.0;
  }

  return state;
}

template State<double, Electron<uint16_t>>
product_state<double, uint16_t>(std::vector<std::string> const &local_states,
                                Electron<uint16_t> const &block);
template State<double, Electron<uint32_t>>
product_state<double, uint32_t>(std::vector<std::string> const &local_states,
                                Electron<uint32_t> const &block);
template State<double, Electron<uint64_t>>
product_state<double, uint64_t>(std::vector<std::string> const &local_states,
                                Electron<uint64_t> const &block);

template State<complex, Electron<uint16_t>>
product_state<complex, uint16_t>(std::vector<std::string> const &local_states,
                                 Electron<uint16_t> const &block);
template State<complex, Electron<uint32_t>>
product_state<complex, uint32_t>(std::vector<std::string> const &local_states,
                                 Electron<uint32_t> const &block);
template State<complex, Electron<uint64_t>>
product_state<complex, uint64_t>(std::vector<std::string> const &local_states,
                                 Electron<uint64_t> const &block);

template <typename coeff_t = complex, typename bit_t = std_bit_t>
State<coeff_t, tJ<bit_t>>
product_state(std::vector<std::string> const &local_states,
              tJ<bit_t> const &block) {
  int n_sites = block.n_sites();
  auto state = zero_state(block);

  if (local_states.size() != n_sites) {
    Log.err("Error creating product state: length of local_states does not "
            "match number of sites in block.");
  }

  if (block.symmetric()) {
    Log.err("Error creating product state: cannot create product state on "
            "symmetric block");
  }

  // Create bit pattern
  bit_t up_pattern = 0;
  bit_t dn_pattern = 0;
  for (int s = 0; s < n_sites; ++s) {
    if (local_states[s] == "Up") {
      up_pattern |= ((bit_t)1 << s);
    } else if (local_states[s] == "Dn") {
      dn_pattern |= ((bit_t)1 << s);
    } else if (local_states[s] == "UpDn") {
      Log.err("Error creating product state: doubly occupied sites not allowed "
              "for t-J block");
    }
    if (local_states[s] != "Emp") {
      Log.err(
          "Error creating product state: invalid local state encountered: {}",
          local_states[s]);
    }
  }

  // using namespace indexing::tj;
  if (block.charge_conserved()) {
    // auto const &indexing = std::get<IndexingNp<bit_t>>(block.indexing());
    // idx_t idx = indexing.index(up_pattern, dn_pattern);
    // state.vector()(idx) = 1.0;
  } else {
    Log.err("Error: tJ block assumes charge conservation");
  }
  return state;
}

template State<double, tJ<uint16_t>>
product_state<double, uint16_t>(std::vector<std::string> const &local_states,
                                tJ<uint16_t> const &block);
template State<double, tJ<uint32_t>>
product_state<double, uint32_t>(std::vector<std::string> const &local_states,
                                tJ<uint32_t> const &block);
template State<double, tJ<uint64_t>>
product_state<double, uint64_t>(std::vector<std::string> const &local_states,
                                tJ<uint64_t> const &block);

template State<complex, tJ<uint16_t>>
product_state<complex, uint16_t>(std::vector<std::string> const &local_states,
                                 tJ<uint16_t> const &block);
template State<complex, tJ<uint32_t>>
product_state<complex, uint32_t>(std::vector<std::string> const &local_states,
                                 tJ<uint32_t> const &block);
template State<complex, tJ<uint64_t>>
product_state<complex, uint64_t>(std::vector<std::string> const &local_states,
                                 tJ<uint64_t> const &block);

} // namespace hydra
