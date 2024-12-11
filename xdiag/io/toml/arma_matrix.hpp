#pragma once

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/extern/toml++/toml.hpp>

namespace xdiag::io {

template <typename T> arma::Mat<T> arma_matrix(toml::node const &node);
template <typename T> toml::array toml_array(arma::Mat<T> const &mat);

} // namespace xdiag::io
