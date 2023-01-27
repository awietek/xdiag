#pragma once

#include <array>
#include <extern/armadillo/armadillo>
#include <extern/toml++/toml.h>
#include <hydra/operators/bond.h>
#include <vector>

namespace hydra::io {

// Read access
using toml_node_view = toml::node_view<toml::node>;
using toml_const_node_view = toml::node_view<const toml::node>;

template <typename T> T get_toml_value(toml_const_node_view const &node);
template <typename T> T get_toml_value(toml_node_view const &node);
template <typename T> T get_toml_value(toml::node const &node);
toml::array get_toml_array(toml_const_node_view const &node);
toml::array get_toml_array(toml_node_view const &node);
toml::array get_toml_array(toml::node const &node);

// toml -> hydra
template <typename T, std::size_t size>
std::array<T, size> toml_array_to_std_array(toml::array const &toml_array);
  
template <typename T>
std::vector<T> toml_array_to_std_vector(toml::array const &toml_array);

template <typename T>
arma::Col<T> toml_array_to_arma_vector(toml::array const &toml_array);

template <typename T>
arma::Mat<T> toml_array_to_arma_matrix(toml::array const &toml_array);

Bond toml_array_to_bond(toml::array const &array);
  
// hydra -> toml
template <typename T>
toml::array std_vector_to_toml_array(std::vector<T> const &value);

template <typename T>
toml::array arma_vector_to_toml_array(arma::Col<T> const &value);

template <typename T>
toml::array arma_matrix_to_toml_array(arma::Mat<T> const &value);

toml::array bond_to_toml_array(Bond const &bond);

} // namespace hydra::io
