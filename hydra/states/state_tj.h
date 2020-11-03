#ifndef HYDRA_STATES_STATE_TJ_H_
#define HYDRA_STATES_STATE_TJ_H_

#include <string>

#include <hydra/common.h>
#include <hydra/qns/qn_tj.h>
#include <hydra/utils/bitops.h>

#include <hydra/combinatorics/up_down_hole.h>
#include <hydra/states/state_spinhalf.h>

#include <iostream>

namespace hydra {

// state for coding tJ models
template <class bit_t = std_bit_t> struct state_tj {
  bit_t ups;
  bit_t dns;
};

// raw index in 3**N Hilbertspace
template <class bit_t = std_bit_t>
inline bit_t rawidx(state_tj<bit_t> state, number_t const &n_sites = 0) {
  bit_t idx = 0;
  bit_t base = 1;
  while (state.ups | state.dns) {
    bit_t local = ((state.ups & (bit_t)1) << 1) | (state.dns & (bit_t)1);
    idx += local * base;
    state.ups >>= 1;
    state.dns >>= 1;
    base *= 3;
  }
  return idx;
}

template <typename bit_t>
std::string String(state_tj<bit_t> state, number_t n_sites);

// Comparison operators
template <class bit_t>
inline bool operator==(state_tj<bit_t> const &s1, state_tj<bit_t> const &s2) {
  return (s1.ups == s2.ups) && (s1.dns == s2.dns);
}
template <class bit_t>
inline bool operator<(state_tj<bit_t> const &s1, state_tj<bit_t> const &s2) {
  return (s1.ups < s2.ups) || ((s1.ups == s2.ups) && (s1.dns < s2.dns));
}
template <class bit_t>
inline bool operator!=(state_tj<bit_t> const &s1, state_tj<bit_t> const &s2) {
  return !(s1 == s2);
}
template <class bit_t>
inline bool operator<=(state_tj<bit_t> const &s1, state_tj<bit_t> const &s2) {
  return s1 < s2 || s1 == s2;
}
template <class bit_t>
inline bool operator>(state_tj<bit_t> const &s1, state_tj<bit_t> const &s2) {
  return !(s1 <= s2);
}
template <class bit_t>
inline bool operator>=(state_tj<bit_t> const &s1, state_tj<bit_t> const &s2) {
  return !(s1 < s2);
}

// get quantum number
template <class bit_t> inline qn_tj QN(state_tj<bit_t> const &s) {
  return qn_tj({utils::popcnt(s.ups), utils::popcnt(s.dns)});
}

// Shift operators
  template <class bit_t, class int_t>
inline state_tj<bit_t> operator<<(state_tj<bit_t> const &s, int_t const &L) {
  return state_tj<bit_t>({s.ups << L, s.dns << L});
}
template <class bit_t, class int_t>
inline state_tj<bit_t> operator>>(state_tj<bit_t> const &s, int_t const &R) {
  return state_tj<bit_t>({s.ups >> R, s.dns >> R});
}
template <class bit_t, class int_t>
inline state_tj<bit_t> &operator<<=(state_tj<bit_t> &s, int_t const &L) {
  s.ups <<= L;
  s.dns <<= L;
  return s;
}
template <class bit_t, class int_t>
inline state_tj<bit_t> &operator>>=(state_tj<bit_t> &s, int_t const &R) {
  s.ups >>= R;
  s.dns >>= R;
  return s;
}

// Substate operators
template <class bit_t>
inline void set_siteval(state_tj<bit_t> &s, number_t const &n,
                        state_tj<bit_t> const &b) {
  utils::sbit(s.ups, n, b.ups);
  utils::sbit(s.dns, n, b.dns);
}

template <class bit_t>
inline state_tj<bit_t> siteval(state_tj<bit_t> const &s, number_t const &n) {
  return state_tj<bit_t>({utils::gbit(s.ups, n), utils::gbit(s.dns, n)});
}
template <class bit_t>
inline state_tj<bit_t> sitevals(state_tj<bit_t> const &s, number_t const &m,
                                number_t const &n) {
  return state_tj<bit_t>(
      {utils::gbits(s.ups, m, n), utils::gbits(s.dns, m, n)});
}

// Combination operator
template <class bit_t>
inline state_tj<bit_t> operator|(state_tj<bit_t> const &s1,
                                 state_tj<bit_t> const &s2) {
  return state_tj<bit_t>({s1.ups | s2.ups, s1.dns | s2.dns});
}

template <class bit_t>
inline state_tj<bit_t> &operator|=(state_tj<bit_t> &s1,
                                   state_tj<bit_t> const &s2) {
  s1.ups |= s2.ups;
  s1.dns |= s2.dns;
  return s1;
}

// Functions to convert spin + hole configuration to up/dn configuration
template <class bit_t = std_bit_t>
inline state_spinhalf<bit_t> down_hole_to_up(state_spinhalf<bit_t> downspins,
                                             state_spinhalf<bit_t> holes) {
  return {combinatorics::down_hole_to_up(downspins.spins, holes.spins)};
}
template <class bit_t = std_bit_t>
inline state_spinhalf<bit_t> up_hole_to_down(state_spinhalf<bit_t> upspins,
                                             state_spinhalf<bit_t> holes) {
  return {combinatorics::up_hole_to_down(upspins.spins, holes.spins)};
}
template <class bit_t = std_bit_t>
inline state_spinhalf<bit_t> up_down_to_hole(state_spinhalf<bit_t> upspins,
                                             state_spinhalf<bit_t> downspins) {
  return {combinatorics::up_down_to_hole(upspins.spins, downspins.spins)};
}
template <class bit_t = std_bit_t>
inline state_spinhalf<bit_t> down_up_to_hole(state_spinhalf<bit_t> downspins,
                                             state_spinhalf<bit_t> upspins) {
  return {combinatorics::down_up_to_hole(downspins.spins, upspins.spins)};
}

} // namespace hydra

#endif
