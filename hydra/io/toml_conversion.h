#pragma once

#include <array>
#include <extern/armadillo/armadillo>
#include <extern/toml++/toml.h>
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

template <typename T, std::size_t size>
std::array<T, size> toml_array_to_std_array(toml::array const &toml_array);

template <typename T>
std::vector<T> toml_array_to_std_vector(toml::array const &toml_array);

template <typename T>
arma::Col<T> toml_array_to_arma_vector(toml::array const &toml_array);

template <typename T>
arma::Mat<T> toml_array_to_arma_matrix(toml::array const &toml_array);

// Write access
template <typename T>
void insert(std::string const &key, T const &value, toml::table &table) {
  table.insert_or_assign(key, value);
}

template <typename T>
void insert_std_vector(std::string const &key, std::vector<T> const &value,
                       toml::table &table);

template <typename T>
void insert_arma_vector(std::string const &key, arma::Col<T> const &value,
                        toml::table &table);

template <typename T>
void insert_arma_matrix(std::string const &key, arma::Mat<T> const &value,
                        toml::table &table);
} // namespace hydra::io
