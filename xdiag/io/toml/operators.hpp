#pragma once

#include <xdiag/extern/toml++/toml.hpp>
#include <xdiag/operators/scalar.hpp>
#include <xdiag/operators/matrix.hpp>
#include <xdiag/operators/coupling.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::io {

Scalar scalar(toml::node const &node);
Coupling coupling(toml::node const &node);
Matrix matrix(toml::node const &node);

Op op(toml::node const &node);
OpSum opsum(toml::node const &node);

} // namespace xdiag::io
