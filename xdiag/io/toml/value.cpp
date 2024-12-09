#include "value.hpp"

#include <xdiag/common.hpp>
#include <xdiag/utils/type_string.hpp>

namespace xdiag::io::toml {

template <typename T> T value(toml::node const &node) try {
  auto val = node.value<T>();
  if (val) {
    return T(*val);
  } else {
    XDIAG_THROW("TOML node cannot be converted to type \"{}\"",
                type_string<T>());
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> complex value<complex>(toml::node const &node) try {
  auto val = node.value<double>();
  auto arr = node.as_array();
  if (val) {
    return complex(*val);
  } else if (arr) {
    auto real_imag = toml_array_to_std_array<double, 2>(*arr);
    return complex{real_imag[0], real_imag[1]};
  } else {
    XDIAG_THROW("TOML node cannot be converted to type \"{}\"",
                type_string<complex>());
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template bool value<bool>(toml::node const &);
template int8_t value<int8_t>(toml::node const &);
template int16_t value<int16_t>(toml::node const &);
template int32_t value<int32_t>(toml::node const &);
template int64_t value<int64_t>(toml::node const &);
template uint8_t value<uint8_t>(toml::node const &);
template uint16_t value<uint16_t>(toml::node const &);
template uint32_t value<uint32_t>(toml::node const &);
template uint64_t value<uint64_t>(toml::node const &);
template double value<double>(toml::node const &);
template std::string value<std::string>(toml::node const &);

} // namespace xdiag::io::toml
