#pragma once

#include <xdiag/extern/toml++/toml.hpp>

namespace xdiag::io::toml {

template <typename T> T value(toml::node const &node);

} // namespace xdiag::io::toml
