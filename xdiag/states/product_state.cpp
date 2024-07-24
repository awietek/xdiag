#include "product_state.hpp"

namespace xdiag {

ProductState::ProductState(std::vector<std::string> const &local_states)
    : local_states_(local_states) {}

std::string const &ProductState::operator[](int64_t i) const {
  return local_states_[i];
}
std::string &ProductState::operator[](int64_t i) { return local_states_[i]; }

int64_t ProductState::n_sites() const { return local_states_.size(); }
void ProductState::operator<<(std::string l) { local_states_.push_back(l); }
ProductState::iterator_t ProductState::begin() const {
  return local_states_.begin();
}
ProductState::iterator_t ProductState::end() const {
  return local_states_.end();
}

State product_state(Block const &block,
                    std::vector<std::string> const &local_state, bool real) {
  return std::visit(
      [&](auto &&block) { return product_state(block, local_state, real); },
      block);
}

template <typename block_t>
State product_state(block_t const &block,
                    std::vector<std::string> const &local_state, bool real) {
  auto state = State(block, real);
  auto pstate = ProductState(local_state);
  fill(state, pstate);
  return state;
}

template State product_state(Spinhalf const &block,
                             std::vector<std::string> const &local_state,
                             bool real);
template State product_state(tJ const &block,
                             std::vector<std::string> const &local_state,
                             bool real);
template State product_state(Electron const &block,
                             std::vector<std::string> const &local_state,
                             bool real);

void fill(State &state, ProductState const &pstate, int64_t col) try {
  if (state.n_sites() != pstate.n_sites()) {
    XDIAG_THROW("State and ProductState do not have the same number of sites");
  } else if (col >= state.n_cols()) {
    XDIAG_THROW("Column index larger than number of columns in State");
  }
  auto const &block = state.block();

  if (state.isreal()) {
    arma::vec v = state.vector(col, false);
    std::visit(
        overload{[&](Spinhalf const &block) { fill(block, v, pstate); },
                 [&](tJ const &block) { fill(block, v, pstate); },
                 [&](Electron const &block) { fill(block, v, pstate); },
#ifdef XDIAG_USE_MPI
                 [&](tJDistributed const &block) { fill(block, v, pstate); },
#endif
                 [&](auto &&) {
                   XDIAG_THROW(
                       "Cannot fill product state for given block (maybe not "
                       "implemented)");
                 }},
        block);
  } else {
    arma::cx_vec v = state.vectorC(col, false);
    std::visit(
        overload{[&](Spinhalf const &block) { fill(block, v, pstate); },
                 [&](tJ const &block) { fill(block, v, pstate); },
                 [&](Electron const &block) { fill(block, v, pstate); },
#ifdef XDIAG_USE_MPI
                 [&](tJDistributed const &block) { fill(block, v, pstate); },
#endif
                 [&](auto &&) {
                   XDIAG_THROW(
                       "Cannot fill product state for given block (maybe not "
                       "implemented)");
                 }},
        block);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename bit_t>
static bit_t spinhalf_bits(ProductState const &pstate, int64_t nup = -1) try {
  bit_t spins = 0;
  int64_t pnup = 0;
  for (int64_t s = 0; s < pstate.n_sites(); ++s) {
    if (pstate[s] == "Up") {
      spins |= ((bit_t)1 << s);
      ++pnup;
    } else {
      if (pstate[s] != "Dn") {
        XDIAG_THROW("Invalid local state encountered");
      }
    }
  }
  if ((nup != pnup) && (nup != -1)) {
    XDIAG_THROW("ProductState does not have correct nup for Spinhalf block");
  }
  return spins;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return 0;
}

template <typename basis_t>
static int64_t spinhalf_index(basis_t const &basis,
                              ProductState const &pstate) try {
  using bit_t = typename basis_t::bit_t;
  bit_t spins = 0;
  if constexpr (basis_t::sz_conserved()) {
    spins = spinhalf_bits<bit_t>(pstate, basis.n_up());
  } else {
    spins = spinhalf_bits<bit_t>(pstate);
  }
  return basis.index(spins);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return 0;
}

template <typename coeff_t>
void fill(Spinhalf const &block, arma::Col<coeff_t> &vec,
          ProductState const &p) try {
  using namespace basis::spinhalf;
  auto const &basis = block.basis();
  int64_t idx = std::visit(
      overload{
          [&](BasisSz<uint32_t> const &b) { return spinhalf_index(b, p); },
          [&](BasisSz<uint64_t> const &b) { return spinhalf_index(b, p); },
          [&](BasisNoSz<uint32_t> const &b) { return spinhalf_index(b, p); },
          [&](BasisNoSz<uint64_t> const &b) { return spinhalf_index(b, p); },
          [&](auto &&) {
            XDIAG_THROW(
                "Cannot create a ProductState for the given type of Basis");
            return (int64_t)-1;
          }},
      basis);
  vec.zeros();
  vec(idx) = 1.0;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void fill(Spinhalf const &, arma::vec &, ProductState const &);
template void fill(Spinhalf const &, arma::cx_vec &, ProductState const &);

template <typename bit_t>
std::pair<bit_t, bit_t> tj_bits(ProductState const &pstate, int64_t nup = -1,
                                int64_t ndn = -1) try {

  bit_t ups = 0;
  bit_t dns = 0;
  int64_t pnup = 0;
  int64_t pndn = 0;
  for (int64_t s = 0; s < pstate.n_sites(); ++s) {
    if (pstate[s] == "Up") {
      ups |= ((bit_t)1 << s);
      ++pnup;
    } else if (pstate[s] == "Dn") {
      dns |= ((bit_t)1 << s);
      ++pndn;
    } else if (pstate[s] == "UpDn") {
      XDIAG_THROW("doubly occupied sites not allowed for t-J block");
      ++pnup;
      ++pndn;
    } else {
      if (pstate[s] != "Emp") {
        XDIAG_THROW("Invalid local state encountered");
      }
    }
  }

  if (((nup != pnup) || (ndn != pndn)) && (nup != -1) && (ndn != -1)) {
    XDIAG_THROW("ProductState does not have correct (nup, ndn) for tJ block");
  }
  return {ups, dns};
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return {0, 0};
}

template <typename basis_t>
static int64_t tj_index(basis_t const &b, ProductState const &pstate) try {
  using bit_t = typename basis_t::bit_t;
  if constexpr (basis_t::np_conserved()) {
    auto [ups, dns] = tj_bits<bit_t>(pstate, b.n_up(), b.n_dn());
    return b.index(ups, dns);
  } else {
    auto [ups, dns] = tj_bits<bit_t>(pstate, -1, -1);
    return b.index(ups, dns);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return 0;
}

template <typename coeff_t>
void fill(tJ const &block, arma::Col<coeff_t> &vec, ProductState const &p) try {
  using namespace basis::tj;
  auto const &basis = block.basis();
  int64_t idx = std::visit(
      overload{
          [&](BasisNp<uint32_t> const &b) { return tj_index(b, p); },
          [&](BasisNp<uint64_t> const &b) { return tj_index(b, p); },
          [&](auto &&) {
            XDIAG_THROW(
                "Cannot create a ProductState for the given type of Basis");
            return (int64_t)-1;
          }},
      basis);
  vec.zeros();
  vec(idx) = 1.0;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void fill(tJ const &, arma::vec &, ProductState const &);
template void fill(tJ const &, arma::cx_vec &, ProductState const &);

template <typename bit_t>
std::pair<bit_t, bit_t> electron_bits(ProductState const &pstate,
                                      int64_t nup = -1, int64_t ndn = -1) try {
  bit_t ups = 0;
  bit_t dns = 0;
  int64_t pnup = 0;
  int64_t pndn = 0;

  for (int64_t s = 0; s < pstate.n_sites(); ++s) {
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
        XDIAG_THROW("Invalid local state encountered");
      }
    }
  }

  if (((nup != pnup) || (ndn != pndn)) && (nup != -1) && (ndn != -1)) {
    XDIAG_THROW(
        "ProductState does not have correct (nup, ndn) for Electron block");
  }
  return {ups, dns};
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return {0, 0};
}

template <typename basis_t>
static int64_t electron_index(basis_t const &b,
                              ProductState const &pstate) try {
  using bit_t = typename basis_t::bit_t;
  if constexpr (basis_t::np_conserved()) {
    auto [ups, dns] = electron_bits<bit_t>(pstate, b.n_up(), b.n_dn());
    return b.index(ups, dns);
  } else {
    auto [ups, dns] = electron_bits<bit_t>(pstate, -1, -1);
    return b.index(ups, dns);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return 0;
}

template <typename coeff_t>
void fill(Electron const &block, arma::Col<coeff_t> &vec,
          ProductState const &p) try {
  using namespace basis::electron;
  auto const &basis = block.basis();
  int64_t idx = std::visit(
      overload{
          [&](BasisNp<uint32_t> const &b) { return electron_index(b, p); },
          [&](BasisNp<uint64_t> const &b) { return electron_index(b, p); },
          [&](BasisNoNp<uint32_t> const &b) { return electron_index(b, p); },
          [&](BasisNoNp<uint64_t> const &b) { return electron_index(b, p); },
          [&](auto &&) {
            XDIAG_THROW(
                "Cannot create a ProductState for the given type of Basis");
            return (int64_t)-1;
          }},
      basis);
  vec.zeros();
  vec(idx) = 1.0;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void fill(Electron const &, arma::vec &, ProductState const &);
template void fill(Electron const &, arma::cx_vec &, ProductState const &);

#ifdef XDIAG_USE_MPI
template <typename basis_t>
static int64_t tj_distributed_index(basis_t const &b,
                                    ProductState const &pstate) try {
  using bit_t = typename basis_t::bit_t;
  if constexpr (basis_t::np_conserved()) {
    auto [ups, dns] = tj_bits<bit_t>(pstate, b.n_up(), b.n_dn());
    return b.index(ups, dns);
  } else {
    auto [ups, dns] = tj_bits<bit_t>(pstate, -1, -1);
    return b.index(ups, dns);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return 0;
}

template <typename coeff_t>
void fill(tJDistributed const &block, arma::Col<coeff_t> &vec,
          ProductState const &p) try {
  using namespace basis::tj_distributed;
  auto const &basis = block.basis();

  // returns invalid_index if state is not on my process
  int64_t idx = std::visit(
      overload{
          [&](BasisNp<uint32_t> const &b) { return tj_index(b, p); },
          [&](BasisNp<uint64_t> const &b) { return tj_index(b, p); },
          [&](auto &&) {
            XDIAG_THROW(
                "Cannot create a ProductState for the given type of Basis");
            return (int64_t)-1;
          }},
      basis);
  vec.zeros();
  if (idx != invalid_index) {
    vec(idx) = 1.0;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void fill(tJDistributed const &, arma::vec &, ProductState const &);
template void fill(tJDistributed const &, arma::cx_vec &, ProductState const &);

#endif

std::ostream &operator<<(std::ostream &out, ProductState const &state) {
  for (auto s : state) {
    out << s;
  }
  out << "\n";
  return out;
}
std::string to_string(ProductState const &state) {
  return to_string_generic(state);
}

} // namespace xdiag
