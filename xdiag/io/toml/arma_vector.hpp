#pragma once

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/extern/toml++/toml.hpp>

namespace xdiag::io {

template <typename T> arma::Col<T> arma_vector(toml::node const &node);
template <typename T> toml::array toml_array(arma::Col<T> const &value);

} // namespace xdiag::io
