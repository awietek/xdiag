#ifndef HYDRA_STATES_STATE_SPINHALF_H_
#define HYDRA_STATES_STATE_SPINHALF_H_

#include <string>

#include <hydra/common.h>
#include <hydra/qns/qn_spinhalf.h>
#include <hydra/utils/bitops.h>

namespace hydra {

// state for coding spinhalf models
template <class bit_t = std_bit_t> struct state_spinhalf {
  bit_t spins;
  explicit operator bit_t() const { return spins; }
};

// raw index in 2**N Hilbertspace
template <class bit_t = std_bit_t>
inline bit_t rawidx(state_spinhalf<bit_t> const &state,
                    number_t const &n_sites = 0) {
  return state.spins;
}

template <typename bit_t>
std::string String(state_spinhalf<bit_t> state, number_t n_sites);

// Comparison operators
template <class bit_t>
inline bool operator==(state_spinhalf<bit_t> const &s1,
                       state_spinhalf<bit_t> const &s2) {
  return s1.spins == s2.spins;
}

template <class bit_t>
inline bool operator<(state_spinhalf<bit_t> const &s1,
                      state_spinhalf<bit_t> const &s2) {
  return s1.spins < s2.spins;
}

template <class bit_t>
inline bool operator!=(state_spinhalf<bit_t> const &s1,
                       state_spinhalf<bit_t> const &s2) {
  return !(s1 == s2);
}

template <class bit_t>
inline bool operator<=(state_spinhalf<bit_t> const &s1,
                       state_spinhalf<bit_t> const &s2) {
  return s1 < s2 || s1 == s2;
}

template <class bit_t>
inline bool operator>(state_spinhalf<bit_t> const &s1,
                      state_spinhalf<bit_t> const &s2) {
  return !(s1 <= s2);
}

template <class bit_t>
inline bool operator>=(state_spinhalf<bit_t> const &s1,
                       state_spinhalf<bit_t> const &s2) {
  return !(s1 < s2);
}

// get quantum number
template <class bit_t> inline qn_spinhalf QN(state_spinhalf<bit_t> const &s) {
  return qn_spinhalf({utils::popcnt(s.spins)});
}

// Shift operators
template <class bit_t, class int_t>
inline state_spinhalf<bit_t> operator<<(state_spinhalf<bit_t> const &s,
                                        int_t const &L) {
  return state_spinhalf<bit_t>{(bit_t)(s.spins << L)};
}
template <class bit_t, class int_t>
inline state_spinhalf<bit_t> operator>>(state_spinhalf<bit_t> const &s,
                                        int_t const &R) {
  return state_spinhalf<bit_t>{(bit_t)(s.spins >> (bit_t)R)};
}
template <class bit_t, class int_t>
inline state_spinhalf<bit_t> &operator<<=(state_spinhalf<bit_t> &s,
                                          int_t const &L) {
  s.spins <<= L;
  return s;
}
template <class bit_t, class int_t>
inline state_spinhalf<bit_t> &operator>>=(state_spinhalf<bit_t> &s,
                                          int_t const &R) {
  s.spins >>= R;
  return s;
}

// Substate operators
template <class bit_t>
inline void set_siteval(state_spinhalf<bit_t> &s, number_t const &n,
                        state_spinhalf<bit_t> const &b) {
  utils::sbit(s.spins, n, b.spins);
}

template <class bit_t>
inline state_spinhalf<bit_t> siteval(state_spinhalf<bit_t> const &s,
                                     number_t const &n) {
  return state_spinhalf<bit_t>({utils::gbit(s.spins, n)});
}
template <class bit_t>
inline state_spinhalf<bit_t> sitevals(state_spinhalf<bit_t> const &s,
                                      number_t const &m, number_t const &n) {
  return state_spinhalf<bit_t>({utils::gbits(s.spins, m, n)});
}

// Combination operator
template <class bit_t>
inline state_spinhalf<bit_t> operator|(state_spinhalf<bit_t> const &s1,
                                       state_spinhalf<bit_t> const &s2) {
  return state_spinhalf<bit_t>({(bit_t)(s1.spins | s2.spins)});
}

template <class bit_t>
inline state_spinhalf<bit_t> &operator|=(state_spinhalf<bit_t> &s1,
                                         state_spinhalf<bit_t> const &s2) {
  s1.spins |= s2.spins;
  return s1;
}

} // namespace hydra

#endif
