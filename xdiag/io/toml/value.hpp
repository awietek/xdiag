// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/extern/toml++/toml.hpp>

namespace xdiag::io {

template <typename T> T value(toml::node const &node);

} // namespace xdiag::io::toml
