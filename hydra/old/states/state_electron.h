#ifndef HYDRA_STATES_STATE_ELECTRON_H_
#define HYDRA_STATES_STATE_ELECTRON_H_

#include <hydra/common.h>
#include <hydra/qns/qn_electron.h>
#include <hydra/utils/bitops.h>

namespace hydra {

// state for coding electron models
template <class bit_t = std_bit_t> struct state_electron {
  bit_t ups;
  bit_t dns;
};

template <typename bit_t>
std::string String(state_electron<bit_t> state, number_t n_sites);

// raw index in 4**N Hilbertspace
template <class bit_t = std_bit_t>
inline bit_t rawidx(state_electron<bit_t> const &state,
                    number_t const &n_sites) {
  return state.ups << n_sites | state.dns;
}

// Comparison operators
template <class bit_t>
inline bool operator==(state_electron<bit_t> const &s1,
                       state_electron<bit_t> const &s2) {
  return (s1.ups == s2.ups) && (s1.dns == s2.dns);
}
template <class bit_t>
inline bool operator<(state_electron<bit_t> const &s1,
                      state_electron<bit_t> const &s2) {
  return (s1.ups < s2.ups) || ((s1.ups == s2.ups) && (s1.dns < s2.dns));
}
template <class bit_t>
inline bool operator!=(state_electron<bit_t> const &s1,
                       state_electron<bit_t> const &s2) {
  return !(s1 == s2);
}
template <class bit_t>
inline bool operator<=(state_electron<bit_t> const &s1,
                       state_electron<bit_t> const &s2) {
  return s1 < s2 || s1 == s2;
}
template <class bit_t>
inline bool operator>(state_electron<bit_t> const &s1,
                      state_electron<bit_t> const &s2) {
  return !(s1 <= s2);
}
template <class bit_t>
inline bool operator>=(state_electron<bit_t> const &s1,
                       state_electron<bit_t> const &s2) {
  return !(s1 < s2);
}

// get quantum number
template <class bit_t> inline qn_electron QN(state_electron<bit_t> const &s) {
  return qn_electron({utils::popcnt(s.ups), utils::popcnt(s.dns)});
}

// Shift operators
  template <class bit_t, class int_t>
inline state_electron<bit_t> operator<<(state_electron<bit_t> const &s,
                                        int_t const &L) {
  return state_electron<bit_t>({s.ups << L, s.dns << L});
}
template <class bit_t, class int_t>
inline state_electron<bit_t> operator>>(state_electron<bit_t> const &s,
                                        int_t const &R) {
  return state_electron<bit_t>({s.ups >> R, s.dns >> R});
}
template <class bit_t, class int_t>
inline state_electron<bit_t> &operator<<=(state_electron<bit_t> &s,
                                          int_t const &L) {
  s.ups <<= L;
  s.dns <<= L;
  return s;
}
template <class bit_t, class int_t>
inline state_electron<bit_t> &operator>>=(state_electron<bit_t> &s,
                                          int_t const &R) {
  s.ups >>= R;
  s.dns >>= R;
  return s;
}

// Substate operators
template <class bit_t>
inline void set_siteval(state_electron<bit_t> &s, number_t const &n,
                        state_electron<bit_t> const &b) {
  utils::sbit(s.ups, n, b.ups);
  utils::sbit(s.dns, n, b.dns);
}
  
template <class bit_t>
inline state_electron<bit_t> siteval(state_electron<bit_t> const &s,
                                     number_t const &n) {
  return state_electron<bit_t>({utils::gbit(s.ups, n), utils::gbit(s.dns, n)});
}
template <class bit_t>
inline state_electron<bit_t> sitevals(state_electron<bit_t> const &s,
                                      number_t const &m, number_t const &n) {
  return state_electron<bit_t>(
      {utils::gbits(s.ups, m, n), utils::gbits(s.dns, m, n)});
}

// Combination operator
template <class bit_t>
inline state_electron<bit_t> operator|(state_electron<bit_t> const &s1,
                                       state_electron<bit_t> const &s2) {
  return state_electron<bit_t>({s1.ups | s2.ups, s1.dns | s2.dns});
}

template <class bit_t>
inline state_electron<bit_t> &operator|=(state_electron<bit_t> &s1,
                                         state_electron<bit_t> const &s2) {
  s1.ups |= s2.ups;
  s1.dns |= s2.dns;
  return s1;
}

} // namespace hydra

#endif
