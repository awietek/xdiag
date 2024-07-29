#include "product_state.hpp"

#include <xdiag/extern/fmt/color.hpp>

namespace xdiag {

ProductState::ProductState(int64_t n_sites) : local_states_(n_sites) {}

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

bool ProductState::operator==(ProductState const &rhs) const {
  return local_states_ == rhs.local_states_;
}
bool ProductState::operator!=(ProductState const &rhs) const {
  return !operator==(rhs);
}

std::ostream &operator<<(std::ostream &out, ProductState const &state) {
  for (int64_t i = state.n_sites() - 1; i >= 0; --i) {
    out << state[i] << " ";
  }
  out << "\n";
  return out;
}
std::string to_string(ProductState const &state, std::string format) try {
  if (format == "plain") {
    return to_string_generic(state);
  } else if (format == "fancy") {
    std::stringstream ss;
    for (int64_t i = state.n_sites() - 1; i >= 0; --i) {
      std::string s = state[i];
      if (s == "Up") {
        // const char *s = u8"\u2B61";
        const char *s = u8"\u2191";
        ss << fmt::format(fg(fmt::color::light_blue), s);
      } else if (s == "Dn") {
        // const char *s = u8"\u2B63";
        const char *s = u8"\u2193";
        ss << fmt::format(fg(fmt::color::orange), s);
      } else if (s == "UpDn") {
        // const char *s = u8"\u2B65";
        const char *s = u8"\u2195";
        ss << fmt::format(fg(fmt::color::red), s);
      } else if (s == "Emp") {
        const char *s = u8"\u25CC";
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
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
