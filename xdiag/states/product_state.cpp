// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "product_state.hpp"

#include <xdiag/extern/fmt/color.hpp>

namespace xdiag {

ProductState::ProductState(int64_t nsites) : local_states_(nsites) {}

ProductState::ProductState(std::vector<std::string> const &local_states)
    : local_states_(local_states) {}

std::string const &ProductState::operator[](int64_t i) const {
  return local_states_[i];
}
std::string &ProductState::operator[](int64_t i) { return local_states_[i]; }

void ProductState::push_back(std::string l) { local_states_.push_back(l); }
int64_t ProductState::size() const { return local_states_.size(); }
int64_t ProductState::nsites() const { return local_states_.size(); }

ProductState::iterator_t ProductState::begin() const {
  return local_states_.begin();
}
ProductState::iterator_t ProductState::end() const {
  return local_states_.end();
}

bool ProductState::operator==(ProductState const &rhs) const {
  return local_states_ == rhs.local_states_;
}
bool ProductState::operator!=(ProductState const &rhs) const {
  return !operator==(rhs);
}

int64_t size(ProductState const &p) { return p.size(); }
int64_t nsites(ProductState const &p) { return p.nsites(); }
std::ostream &operator<<(std::ostream &out, ProductState const &state) {
  for (int64_t i = state.size() - 1; i >= 0; --i) {
    out << state[i] << " ";
  }
  return out;
}
std::string to_string(ProductState const &state, std::string format) try {
  if (format == "plain") {
    return to_string_generic(state);
  } else if (format == "fancy") {
    std::stringstream ss;
    for (int64_t i = state.size() - 1; i >= 0; --i) {
      std::string s = state[i];
      if (s == "Up") {
        // const char *s = u8"\u2B61";
        const char *s = "\u2191";
        ss << fmt::format(fg(fmt::color::light_blue), s);
      } else if (s == "Dn") {
        // const char *s = u8"\u2B63";
        const char *s = "\u2193";
        ss << fmt::format(fg(fmt::color::orange), s);
      } else if (s == "UpDn") {
        // const char *s = u8"\u2B65";
        const char *s = "\u2195";
        ss << fmt::format(fg(fmt::color::red), s);
      } else if (s == "Emp") {
        const char *s = "\u25CC";
        ss << fmt::format(fg(fmt::color::gray), s);
      } else {
        XDIAG_THROW(fmt::format("Unknown local state for fancy formatting "
                                "style of ProductState: {}",
                                s));
      }
    }
    return ss.str();
  } else {
    XDIAG_THROW(
        fmt::format("Unknown formatting style for ProductState: {}", format));
  }
}
XDIAG_CATCH

template <typename bit_t>
void to_product_state_spinhalf(bit_t spins, ProductState &pstate) {
  for (int64_t i = 0; i < pstate.size(); ++i) {
    if (spins & 1) {
      pstate[i] = "Up";
    } else {
      pstate[i] = "Dn";
    }
    spins >>= 1;
  }
}
template void to_product_state_spinhalf(uint32_t, ProductState &);
template void to_product_state_spinhalf(uint64_t, ProductState &);

template <typename bit_t>
void to_product_state_tj(bit_t ups, bit_t dns, ProductState &pstate) {
  for (int64_t i = 0; i < pstate.size(); ++i) {
    if (ups & 1) {
      pstate[i] = "Up";
    } else {
      if (dns & 1) {
        pstate[i] = "Dn";
      } else {
        pstate[i] = "Emp";
      }
    }
    ups >>= 1;
    dns >>= 1;
  }
}
template void to_product_state_tj(uint32_t ups, uint32_t dns,
                                  ProductState &pstate);
template void to_product_state_tj(uint64_t ups, uint64_t dns,
                                  ProductState &pstate);

template <typename bit_t>
void to_product_state_electron(bit_t ups, bit_t dns, ProductState &pstate) {
  for (int64_t i = 0; i < pstate.size(); ++i) {
    if (ups & 1) {
      if (dns & 1) {
        pstate[i] = "UpDn";
      } else {
        pstate[i] = "Up";
      }
    } else {
      if (dns & 1) {
        pstate[i] = "Dn";
      } else {
        pstate[i] = "Emp";
      }
    }
    ups >>= 1;
    dns >>= 1;
  }
}
template void to_product_state_electron(uint32_t ups, uint32_t dns,
                                        ProductState &pstate);
template void to_product_state_electron(uint64_t ups, uint64_t dns,
                                        ProductState &pstate);

template <typename bit_t>
bit_t to_bits_spinhalf(ProductState const &pstate) try {
  bit_t spins = 0;
  for (int64_t s = 0; s < pstate.size(); ++s) {
    std::string val = pstate[s];
    if (val == "Up") {
      spins |= ((bit_t)1 << s);
    } else {
      if (val != "Dn") {
        XDIAG_THROW(fmt::format("Invalid local state encountered: {}", val));
      }
    }
  }
  return spins;
}
XDIAG_CATCH

template uint32_t to_bits_spinhalf(ProductState const &);
template uint64_t to_bits_spinhalf(ProductState const &);

template <typename bit_t>
std::pair<bit_t, bit_t> to_bits_tj(ProductState const &pstate) try {
  bit_t ups = 0;
  bit_t dns = 0;
  for (int64_t s = 0; s < pstate.size(); ++s) {
    std::string val = pstate[s];
    if (val == "Up") {
      ups |= ((bit_t)1 << s);
    } else if (val == "Dn") {
      dns |= ((bit_t)1 << s);
    } else {
      if (val != "Emp") {
        XDIAG_THROW(fmt::format("Invalid local state encountered: {}", val));
      }
    }
  }
  return {ups, dns};
}
XDIAG_CATCH

template std::pair<uint32_t, uint32_t> to_bits_tj(ProductState const &);
template std::pair<uint64_t, uint64_t> to_bits_tj(ProductState const &);

template <typename bit_t>
std::pair<bit_t, bit_t> to_bits_electron(ProductState const &pstate) try {
  bit_t ups = 0;
  bit_t dns = 0;
  for (int64_t s = 0; s < pstate.size(); ++s) {
    std::string val = pstate[s];
    if (val == "Up") {
      ups |= ((bit_t)1 << s);
    } else if (val == "Dn") {
      dns |= ((bit_t)1 << s);
    } else if (val == "UpDn") {
      ups |= ((bit_t)1 << s);
      dns |= ((bit_t)1 << s);
    } else {
      if (val != "Emp") {
        XDIAG_THROW(fmt::format("Invalid local state encountered: {}", val));
      }
    }
  }
  return {ups, dns};
}
XDIAG_CATCH

template std::pair<uint32_t, uint32_t> to_bits_electron(ProductState const &);
template std::pair<uint64_t, uint64_t> to_bits_electron(ProductState const &);

} // namespace xdiag
