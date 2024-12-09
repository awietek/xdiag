#include "std_vector.hpp"

#include <xdiag/common.hpp>
#include <xdiag/io/toml/value.hpp>
#include <xdiag/utils/type_string.hpp>

namespace xdiag::io {

template <typename T> std::vector<T> std_vector(toml::node const &node) try {
  auto array = node.as_array();
  if (array) {
    std::size_t size = array->size();
    std::vector<T> vector(size);
    for (int i = 0; i < size; ++i) {
      vector[i] = value<T>(array[i]);
    }
    return vector;
  } else {
    XDIAG_THROW(fmt::format(
        "TOML node cannot be converted to \"{}\". Node is not an array.",
        utils::type_string<std::vector<T>>()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template std::vector<int8_t> std_vector<int8_t>(toml::node const &);
template std::vector<int16_t> std_vector<int16_t>(toml::node const &);
template std::vector<int32_t> std_vector<int32_t>(toml::node const &);
template std::vector<int64_t> std_vector<int64_t>(toml::node const &);
template std::vector<uint8_t> std_vector<uint8_t>(toml::node const &);
template std::vector<uint16_t> std_vector<uint16_t>(toml::node const &);
template std::vector<uint32_t> std_vector<uint32_t>(toml::node const &);
template std::vector<uint64_t> std_vector<uint64_t>(toml::node const &);
template std::vector<double> std_vector<double>(toml::node const &);
template std::vector<complex> std_vector<complex>(toml::node const &);
template std::vector<std::string> std_vector<std::string>(toml::node const &);

} // namespace xdiag::io
