#include "representative.hpp"

#include <limits>

namespace xdiag::symmetries {

template <typename bit_t, typename group_action_t>
bool is_representative(bit_t state, group_action_t const &group_action) {
  for (int64_t sym = 0; sym < group_action.nsymmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < state) {
      return false;
    }
  }
  return true;
}

template <typename bit_t, typename group_action_t>
bit_t representative(bit_t state, group_action_t const &group_action) {
  if (group_action.nsymmetries() == 0) {
    return state;
  }

  bit_t rep = std::numeric_limits<bit_t>::max();
  for (int64_t sym = 0; sym < group_action.nsymmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < rep) {
      rep = tstate;
    }
  }
  return rep;
}

template <typename bit_t, class group_action_t>
bit_t representative_subset(bit_t state, group_action_t const &group_action,
                            gsl::span<int64_t const> syms) {
  if (syms.size() == 0) {
    return state;
  }
  bit_t rep = std::numeric_limits<bit_t>::max();
  for (int64_t sym : syms) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < rep) {
      rep = tstate;
    }
  }
  return rep;
}

template <typename bit_t, class group_action_t>
inline std::pair<bit_t, int64_t>
representative_sym(bit_t state, group_action_t const &group_action) {
  if (group_action.n_symmetries() == 0) {
    return {state, 0};
  }
  bit_t rep = std::numeric_limits<bit_t>::max();
  int64_t rep_sym = 0;
  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < rep) {
      rep = tstate;
      rep_sym = sym;
    }
  }
  return {rep, rep_sym};
}

template <typename bit_t, class group_action_t>
std::pair<bit_t, int64_t>
representative_sym_subset(bit_t state, group_action_t const &group_action,
                          gsl::span<int64_t const> syms) {
  if (syms.size() == 0) {
    return {state, 0};
  }
  bit_t rep = std::numeric_limits<bit_t>::max();
  int64_t rep_sym = 0;
  for (auto sym : syms) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < rep) {
      rep = tstate;
      rep_sym = sym;
    }
  }
  return {rep, rep_sym};
}

template <typename bit_t, class group_action_t>
std::pair<bit_t, std::vector<int64_t>>
representative_syms(bit_t state, group_action_t const &group_action) {
  if (group_action.n_symmetries() == 0) {
    return {state, std::vector<int64_t>()};
  }

  bit_t rep = representative(state, group_action);
  std::vector<int64_t> rep_syms;
  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate == rep) {
      rep_syms.push_back(sym);
    }
  }
  return {rep, rep_syms};
}

} // namespace xdiag::symmetries
