#include "product_state.h"

#include <hydra/indexing/indexing_variants.h>

namespace hydra {

ProductState::ProductState(std::vector<std::string> const &local_states)
    : local_states_(local_states) {}

template <typename coeff_t>
void fill(ProductState const &pstate, State<coeff_t> &state) {
  auto &vector = state.vector();
  std::visit(overloaded{[&pstate, &vector](auto &&block) {
               fill(pstate, block, vector);
             }},
             state.block());
}
template void fill(ProductState const &pstate, State<double> &state);
template void fill(ProductState const &pstate, State<complex> &state);

template <class block_t>
void check_product_state_n_sites(block_t const &block,
                                 ProductState const &pstate) {
  if (block.n_sites() != pstate.n_sites()) {
    Log.err("Error creating product state: number of sites of product state "
            "does not match number of sites in block.");
  }
}

template <class block_t>
void check_product_state_not_symmetric(block_t const &block) {
  if (block.symmetric()) {
    Log.err("Error creating product state: cannot create product state on "
            "a symmetric block");
  }
}

template <typename bit_t>
bit_t get_spinhalf_spins(ProductState const &pstate, int nup = -1) {
  bit_t spins = 0;
  int pnup = 0;
  for (int s = 0; s < pstate.n_sites(); ++s) {
    if (pstate[s] == "Up") {
      spins |= ((bit_t)1 << s);
      ++pnup;
    } else {
      if (pstate[s] != "Dn") {
        Log.err(
            "Error creating product state: invalid local state encountered: {}",
            pstate[s]);
      }
    }
  }

  if ((nup != pnup) && (nup != -1)) {
    Log.err("Error creating product state: number of up spins incompatible "
            "with block");
  }
  return spins;
}

template <typename coeff_t>
void fill(ProductState const &pstate, Spinhalf const &block,
          arma::Col<coeff_t> &vector) {
  using namespace indexing::spinhalf;
  check_product_state_n_sites(block, pstate);
  check_product_state_not_symmetric(block);

  vector.zeros();
  std::visit(overloaded{
                 [&pstate, &vector](IndexingSz<uint16_t> const &indexing) {
                   uint16_t spins =
                       get_spinhalf_spins<uint16_t>(pstate, indexing.n_up());
                   idx_t idx = indexing.index(spins);
                   vector(idx) = 1.0;
                 },
                 [&pstate, &vector](IndexingSz<uint32_t> const &indexing) {
                   uint32_t spins =
                       get_spinhalf_spins<uint32_t>(pstate, indexing.n_up());
                   idx_t idx = indexing.index(spins);
                   vector(idx) = 1.0;
                 },
                 [&pstate, &vector](IndexingSz<uint64_t> const &indexing) {
                   uint64_t spins =
                       get_spinhalf_spins<uint64_t>(pstate, indexing.n_up());
                   idx_t idx = indexing.index(spins);
                   vector(idx) = 1.0;
                 },
                 [&pstate, &vector](IndexingNoSz<uint16_t> const &indexing) {
                   uint16_t spins = get_spinhalf_spins<uint16_t>(pstate);
                   idx_t idx = indexing.index(spins);
                   vector(idx) = 1.0;
                 },
                 [&pstate, &vector](IndexingNoSz<uint32_t> const &indexing) {
                   uint32_t spins = get_spinhalf_spins<uint32_t>(pstate);
                   idx_t idx = indexing.index(spins);
                   vector(idx) = 1.0;
                 },
                 [&pstate, &vector](IndexingNoSz<uint64_t> const &indexing) {
                   uint64_t spins = get_spinhalf_spins<uint16_t>(pstate);
                   idx_t idx = indexing.index(spins);
                   vector(idx) = 1.0;
                 },
                 [](auto const &) {
                   Log.err("Error creating product state: Invalid Indexing");
                 },
             },
             block.indexing());
}

template void fill(ProductState const &pstate, Spinhalf const &block,
                   arma::Col<double> &vector);
template void fill(ProductState const &pstate, Spinhalf const &block,
                   arma::Col<complex> &vector);

template <typename bit_t>
std::pair<bit_t, bit_t> get_tj_spins(ProductState const &pstate, int nup = -1,
                                     int ndn = -1) {

  bit_t ups = 0;
  bit_t dns = 0;
  int pnup = 0;
  int pndn = 0;
  for (int s = 0; s < pstate.n_sites(); ++s) {
    if (pstate[s] == "Up") {
      ups |= ((bit_t)1 << s);
      ++pnup;
    } else if (pstate[s] == "Dn") {
      dns |= ((bit_t)1 << s);
      ++pndn;
    } else if (pstate[s] == "UpDn") {
      Log.err("Error creating product state: doubly occupied sites not allowed "
              "for t-J block");
      ++pnup;
      ++pndn;
    } else {
      if (pstate[s] != "Emp") {
        Log.err(
            "Error creating product state: invalid local state encountered: {}",
            pstate[s]);
      }
    }
  }

  if ((nup != pnup) && (nup != -1) && (ndn != -1)) {
    Log.err("Error creating product state: number of up/dn spins incompatible "
            "with block");
  }

  return {ups, dns};
}

template <typename coeff_t>
void fill(ProductState const &pstate, tJ const &block,
          arma::Col<coeff_t> &vector) {
  using namespace indexing::tj;
  check_product_state_n_sites(block, pstate);
  check_product_state_not_symmetric(block);
  vector.zeros();
  std::visit(overloaded{
                 [&pstate, &vector](IndexingNp<uint16_t> const &indexing) {
                   auto [ups, dns] = get_tj_spins<uint16_t>(
                       pstate, indexing.n_up(), indexing.n_dn());
                   idx_t idx = indexing.index(ups, dns);
                   vector(idx) = 1.0;
                 },
                 [&pstate, &vector](IndexingNp<uint32_t> const &indexing) {
                   auto [ups, dns] = get_tj_spins<uint32_t>(
                       pstate, indexing.n_up(), indexing.n_dn());
                   idx_t idx = indexing.index(ups, dns);
                   vector(idx) = 1.0;
                 },
                 [&pstate, &vector](IndexingNp<uint64_t> const &indexing) {
                   auto [ups, dns] = get_tj_spins<uint64_t>(
                       pstate, indexing.n_up(), indexing.n_dn());
                   idx_t idx = indexing.index(ups, dns);
                   vector(idx) = 1.0;
                 },
                 [](auto const &) {
                   Log.err("Error creating product state: Invalid Indexing");
                 },
             },
             block.indexing());
}

template void fill(ProductState const &pstate, tJ const &block,
                   arma::Col<double> &vector);
template void fill(ProductState const &pstate, tJ const &block,
                   arma::Col<complex> &vector);

template <typename bit_t>
std::pair<bit_t, bit_t> get_electron_spins(ProductState const &pstate,
                                           int nup = -1, int ndn = -1) {

  bit_t ups = 0;
  bit_t dns = 0;
  int pnup = 0;
  int pndn = 0;

  for (int s = 0; s < pstate.n_sites(); ++s) {
    if (pstate[s] == "Up") {
      ups |= ((bit_t)1 << s);
      ++pnup;
    } else if (pstate[s] == "Dn") {
      dns |= ((bit_t)1 << s);
      ++pndn;
    } else if (pstate[s] == "UpDn") {
      ups |= ((bit_t)1 << s);
      dns |= ((bit_t)1 << s);
      ++pnup;
      ++pndn;
    } else {
      if (pstate[s] != "Emp") {
        Log.err(
            "Error creating product state: invalid local state encountered: {}",
            pstate[s]);
      }
    }
  }

  if ((nup != pnup) && (nup != -1) && (ndn != -1)) {
    Log.err("Error creating product state: number of up/dn spins incompatible "
            "with block");
  }

  return {ups, dns};
}

template <typename coeff_t>
void fill(ProductState const &pstate, Electron const &block,
          arma::Col<coeff_t> &vector) {
  using namespace indexing::electron;
  check_product_state_n_sites(block, pstate);
  check_product_state_not_symmetric(block);
  vector.zeros();
  std::visit(overloaded{
                 [&pstate, &vector](IndexingNp<uint16_t> const &indexing) {
                   auto [ups, dns] = get_electron_spins<uint16_t>(
                       pstate, indexing.n_up(), indexing.n_dn());
                   idx_t idx = indexing.index(ups, dns);
                   vector(idx) = 1.0;
                 },
                 [&pstate, &vector](IndexingNp<uint32_t> const &indexing) {
                   auto [ups, dns] = get_electron_spins<uint32_t>(
                       pstate, indexing.n_up(), indexing.n_dn());
                   idx_t idx = indexing.index(ups, dns);
                   vector(idx) = 1.0;
                 },
                 [&pstate, &vector](IndexingNp<uint64_t> const &indexing) {
                   auto [ups, dns] = get_electron_spins<uint64_t>(
                       pstate, indexing.n_up(), indexing.n_dn());
                   idx_t idx = indexing.index(ups, dns);
                   vector(idx) = 1.0;
                 },
                 [&pstate, &vector](IndexingNoNp<uint16_t> const &indexing) {
                   auto [ups, dns] = get_electron_spins<uint16_t>(pstate);
                   idx_t idx = indexing.index(ups, dns);
                   vector(idx) = 1.0;
                 },
                 [&pstate, &vector](IndexingNoNp<uint32_t> const &indexing) {
                   auto [ups, dns] = get_electron_spins<uint32_t>(pstate);
                   idx_t idx = indexing.index(ups, dns);
                   vector(idx) = 1.0;
                 },
                 [&pstate, &vector](IndexingNoNp<uint64_t> const &indexing) {
                   auto [ups, dns] = get_electron_spins<uint64_t>(pstate);
                   idx_t idx = indexing.index(ups, dns);
                   vector(idx) = 1.0;
                 },
                 [](auto const &) {
                   Log.err("Error creating product state: Invalid Indexing");
                 },
             },
             block.indexing());
}

template void fill(ProductState const &pstate, Electron const &block,
                   arma::Col<double> &vector);
template void fill(ProductState const &pstate, Electron const &block,
                   arma::Col<complex> &vector);
} // namespace hydra
