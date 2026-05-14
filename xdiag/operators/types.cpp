// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "types.hpp"

#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag {

bool is_known_type(std::string const &type) {
  return op_registry.find(type) != op_registry.end();
}

bool is_real_type(std::string const &type) {
  auto it = op_registry.find(type);
  return (it != op_registry.end()) && it->second.is_real;
}

int64_t nsites_of_type(std::string const &type) try {
  auto it = op_registry.find(type);
  if (it == op_registry.end()) {
    XDIAG_THROW(fmt::format("Op type \"{}\" is unknown. Known types: {}",
                            type, known_types_string()));
  }
  return it->second.nsites;
}
XDIAG_CATCH

OpTypeInfo const &info_of_type(std::string const &type) try {
  auto it = op_registry.find(type);
  if (it == op_registry.end()) {
    XDIAG_THROW(fmt::format("Op type \"{}\" is unknown. Known types: {}",
                            type, known_types_string()));
  }
  return it->second;
}
XDIAG_CATCH

std::string known_types_string() {
  std::string str;
  for (auto const &[type, info] : op_registry) {
    str += fmt::format("\"{}\", ", type);
  }
  if (!str.empty())
    str.resize(str.size() - 2);
  return str;
}

} // namespace xdiag
