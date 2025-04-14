// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/extern/toml++/toml.hpp>
#include <xdiag/operators/coupling.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/matrix.hpp>
#include <xdiag/utils/scalar.hpp>

namespace xdiag::io {

Scalar scalar(toml::node const &node);
Coupling coupling(toml::node const &node);
Matrix matrix(toml::node const &node);

Op op(toml::node const &node);
OpSum opsum(toml::node const &node);

toml::array toml_array(Op const &op);
toml::array toml_array(OpSum const &ops);
toml::table toml_table(OpSum const &ops);

} // namespace xdiag::io
