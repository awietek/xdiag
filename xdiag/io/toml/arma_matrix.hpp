#pragma once

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/extern/toml++/toml.hpp>

namespace xdiag::io {

template <typename T> arma::Mat<T> arma_matrix(toml::node const &node);

} // namespace xdiag::io
