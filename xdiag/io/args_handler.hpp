#pragma once

#include <map>
#include <string>
#include <variant>
#include <vector>

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/common.hpp>

namespace xdiag::io {

using args_t =
    std::variant<std::string, bool, int8_t, int16_t, int32_t, int64_t, uint8_t,
                 uint16_t, uint32_t, uint64_t, double, complex>;


class ArgsHandler {
public:
  ArgsHandler(std::string key, std::map<std::string, io::args_t> &args);
  ArgsHandler(ArgsHandler const &) = delete;
  ArgsHandler &operator=(ArgsHandler const &) = delete;

  template <class data_t> data_t as() const;
  template <class data_t> data_t as(data_t const &defaultt) const;
  template <class data_t> void operator=(data_t const &data);

private:
  std::string key_;
  std::map<std::string, io::args_t> &args_;
};

} // namespace xdiag::io
