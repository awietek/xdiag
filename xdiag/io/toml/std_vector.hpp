// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>
#include <xdiag/extern/toml++/toml.hpp>

namespace xdiag::io {

template <typename T> std::vector<T> std_vector(toml::node const &node);
template <typename T> toml::array toml_array(std::vector<T> const &vec);

} // namespace xdiag::io
