#include "toml_conversion.h"

#include <hydra/common.h>
#include <hydra/utils/logger.h>

namespace hydra::io {

template <typename T> T get_toml_value(toml_const_node_view const &node) {
  auto val = node.value<T>();
  if (val) {
    return T(*val);
  } else {
    Log.err("Error parsing toml: node does not contain a value!");
    return T();
  }
}

template <> complex get_toml_value<complex>(toml_const_node_view const &node) {
  auto val = node.value<double>();
  auto arr = node.as_array();
  if (val) {
    return complex(*val);
  } else if (arr) {
    auto real_imag = toml_array_to_std_array<double, 2>(*arr);
    return complex{real_imag[0], real_imag[1]};
  } else {
    Log.err("Error parsing toml: failed trying to parse complex number!");
  }
}

template bool get_toml_value<bool>(toml_const_node_view const &);
template int8_t get_toml_value<int8_t>(toml_const_node_view const &);
template int16_t get_toml_value<int16_t>(toml_const_node_view const &);
template int32_t get_toml_value<int32_t>(toml_const_node_view const &);
template int64_t get_toml_value<int64_t>(toml_const_node_view const &);
template uint8_t get_toml_value<uint8_t>(toml_const_node_view const &);
template uint16_t get_toml_value<uint16_t>(toml_const_node_view const &);
template uint32_t get_toml_value<uint32_t>(toml_const_node_view const &);
template uint64_t get_toml_value<uint64_t>(toml_const_node_view const &);
template double get_toml_value<double>(toml_const_node_view const &);
template std::string get_toml_value<std::string>(toml_const_node_view const &);

template <typename T> T get_toml_value(toml_node_view const &node) {
  auto val = node.value<T>();
  if (val) {
    return T(*val);
  } else {
    Log.err("Error parsing toml: node does not contain a value!");
    return T();
  }
}

template <> complex get_toml_value<complex>(toml_node_view const &node) {
  auto val = node.value<double>();
  auto arr = node.as_array();
  if (val) {
    return complex(*val);
  } else if (arr) {
    auto real_imag = toml_array_to_std_array<double, 2>(*arr);
    return complex{real_imag[0], real_imag[1]};
  } else {
    Log.err("Error parsing toml: failed trying to parse complex number!");
  }
}

template bool get_toml_value<bool>(toml_node_view const &);
template int8_t get_toml_value<int8_t>(toml_node_view const &);
template int16_t get_toml_value<int16_t>(toml_node_view const &);
template int32_t get_toml_value<int32_t>(toml_node_view const &);
template int64_t get_toml_value<int64_t>(toml_node_view const &);
template uint8_t get_toml_value<uint8_t>(toml_node_view const &);
template uint16_t get_toml_value<uint16_t>(toml_node_view const &);
template uint32_t get_toml_value<uint32_t>(toml_node_view const &);
template uint64_t get_toml_value<uint64_t>(toml_node_view const &);
template double get_toml_value<double>(toml_node_view const &);
template std::string get_toml_value<std::string>(toml_node_view const &);

template <typename T> T get_toml_value(toml::node const &node) {
  auto val = node.value<T>();
  if (val) {
    return T(*val);
  } else {
    Log.err("Error parsing toml: node does not contain a value!");
    return T();
  }
}

template <> complex get_toml_value<complex>(toml::node const &node) {
  auto val = node.value<double>();
  auto arr = node.as_array();
  if (val) {
    return complex(*val);
  } else if (arr) {
    auto real_imag = toml_array_to_std_array<double, 2>(*arr);
    return complex{real_imag[0], real_imag[1]};
  } else {
    Log.err("Error parsing toml: failed trying to parse complex number!");
  }
}
template bool get_toml_value<bool>(toml::node const &);
template int8_t get_toml_value<int8_t>(toml::node const &);
template int16_t get_toml_value<int16_t>(toml::node const &);
template int32_t get_toml_value<int32_t>(toml::node const &);
template int64_t get_toml_value<int64_t>(toml::node const &);
template uint8_t get_toml_value<uint8_t>(toml::node const &);
template uint16_t get_toml_value<uint16_t>(toml::node const &);
template uint32_t get_toml_value<uint32_t>(toml::node const &);
template uint64_t get_toml_value<uint64_t>(toml::node const &);
template double get_toml_value<double>(toml::node const &);
template complex get_toml_value<complex>(toml::node const &);
template std::string get_toml_value<std::string>(toml::node const &);

toml::array get_toml_array(toml_const_node_view const &node) {
  auto arr = node.as_array();
  if (arr) {
    return *arr;
  } else {
    Log.err("Error parsing toml: node does not contain an array!");
    return toml::array();
  }
}

toml::array get_toml_array(toml_node_view const &node) {
  auto arr = node.as_array();
  if (arr) {
    return *arr;
  } else {
    Log.err("Error parsing toml: node does not contain an array!");
    return toml::array();
  }
}

toml::array get_toml_array(toml::node const &node) {
  auto arr = node.as_array();
  if (arr) {
    return *arr;
  } else {
    Log.err("Error parsing toml: node does not contain an array!");
    return toml::array();
  }
}

template <typename T, std::size_t size>
std::array<T, size> toml_array_to_std_array(toml::array const &toml_array) {
  if (size != toml_array.size()) {
    Log.err("Error reading TOML array to std::array: toml::array size ({}) "
            "exceeds std::array size ({}) !",
            toml_array.size(), size);
  }

  std::array<T, size> array;
  for (std::size_t i = 0; i < size; ++i) {
    array[i] = get_toml_value<T>(toml_array[i]);
  }
  return array;
}

template std::array<double, 2>
toml_array_to_std_array<double, 2>(toml::array const &);

template <typename T>
std::vector<T> toml_array_to_std_vector(toml::array const &toml_array) {
  std::size_t size = toml_array.size();
  std::vector<T> vector(size);
  for (std::size_t i = 0; i < size; ++i) {
    vector[i] = get_toml_value<T>(toml_array[i]);
  }
  return vector;
}

template std::vector<int8_t>
toml_array_to_std_vector<int8_t>(toml::array const &);
template std::vector<int16_t>
toml_array_to_std_vector<int16_t>(toml::array const &);
template std::vector<int32_t>
toml_array_to_std_vector<int32_t>(toml::array const &);
template std::vector<int64_t>
toml_array_to_std_vector<int64_t>(toml::array const &);
template std::vector<uint8_t>
toml_array_to_std_vector<uint8_t>(toml::array const &);
template std::vector<uint16_t>
toml_array_to_std_vector<uint16_t>(toml::array const &);
template std::vector<uint32_t>
toml_array_to_std_vector<uint32_t>(toml::array const &);
template std::vector<uint64_t>
toml_array_to_std_vector<uint64_t>(toml::array const &);
template std::vector<double>
toml_array_to_std_vector<double>(toml::array const &);
template std::vector<complex>
toml_array_to_std_vector<complex>(toml::array const &);
template std::vector<std::string>
toml_array_to_std_vector<std::string>(toml::array const &);

template <typename T>
arma::Col<T> toml_array_to_arma_vector(toml::array const &toml_array) {
  std::size_t size = toml_array.size();
  arma::Col<T> vector(size);
  for (arma::uword i = 0; i < size; ++i) {
    vector[i] = get_toml_value<T>(toml_array[i]);
  }
  return vector;
}

template arma::Col<double>
toml_array_to_arma_vector<double>(toml::array const &);
template arma::Col<complex>
toml_array_to_arma_vector<complex>(toml::array const &);
template arma::Col<arma::sword>
toml_array_to_arma_vector<arma::sword>(toml::array const &);
template arma::Col<arma::uword>
toml_array_to_arma_vector<arma::uword>(toml::array const &);

template <typename T>
arma::Mat<T> toml_array_to_arma_matrix(toml::array const &toml_array) {
  toml::array const &rows = toml_array;
  std::size_t m = rows.size();
  std::size_t n = 0;
  if (m > 0) {
    std::vector<T> vector;

    for (std::size_t i = 0; i < m; ++i) {
      auto row_entries_opt = rows[i].as_array();
      if (row_entries_opt) {
        toml::array row_entries = *row_entries_opt;
        if (i == 0) {
          n = row_entries.size();
          vector.resize(m * n);
        }
        if (row_entries.size() != n) {
          Log.err("Error reading TOML array to arma::Matrix: not all rows "
                  "have same length");
        }

        for (std::size_t j = 0; j < n; ++j) {
          vector[i + m * j] = get_toml_value<T>(row_entries[j]);
        }

      } else {
        Log.err(
            "Error reading TOML array to arma::Matrix: Cannot parse row {}!",
            i);
      }
    } // for (int i = 0; i < m; ++i)
    return arma::Mat<T>(vector.data(), m, n);
  } else {
    return arma::Mat<T>();
  }
}

template arma::Mat<double>
toml_array_to_arma_matrix<double>(toml::array const &);
template arma::Mat<complex>
toml_array_to_arma_matrix<complex>(toml::array const &);
template arma::Mat<arma::sword>
toml_array_to_arma_matrix<arma::sword>(toml::array const &);
template arma::Mat<arma::uword>
toml_array_to_arma_matrix<arma::uword>(toml::array const &);

template <typename T>
void insert_std_vector(std::string const &key, std::vector<T> const &value,
                       toml::table &table) {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(x);
  }
  insert(key, arr, table);
}
template <>
void insert_std_vector(std::string const &key,
                       std::vector<uint64_t> const &value, toml::table &table) {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back((int64_t)x);
  }
  insert(key, arr, table);
}

template <>
void insert_std_vector(std::string const &key,
                       std::vector<complex> const &value, toml::table &table) {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(toml::array{x.real(), x.imag()});
  }
  insert(key, arr, table);
}

template void insert_std_vector(std::string const &,
                                std::vector<int8_t> const &, toml::table &);
template void insert_std_vector(std::string const &,
                                std::vector<int16_t> const &, toml::table &);
template void insert_std_vector(std::string const &,
                                std::vector<int32_t> const &, toml::table &);
template void insert_std_vector(std::string const &,
                                std::vector<int64_t> const &, toml::table &);
template void insert_std_vector(std::string const &,
                                std::vector<uint8_t> const &, toml::table &);
template void insert_std_vector(std::string const &,
                                std::vector<uint16_t> const &, toml::table &);
template void insert_std_vector(std::string const &,
                                std::vector<uint32_t> const &, toml::table &);
template void insert_std_vector(std::string const &,
                                std::vector<double> const &, toml::table &);
template void insert_std_vector(std::string const &,
                                std::vector<std::string> const &,
                                toml::table &);

template <typename T>
void insert_arma_vector(std::string const &key, arma::Col<T> const &value,
                        toml::table &table) {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(x);
  }
  insert(key, arr, table);
}
template <>
void insert_arma_vector(std::string const &key,
                        arma::Col<arma::uword> const &value,
                        toml::table &table) {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back((int64_t)x);
  }
  insert(key, arr, table);
}
template <>
void insert_arma_vector(std::string const &key, arma::Col<complex> const &value,
                        toml::table &table) {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(toml::array{x.real(), x.imag()});
  }
  insert(key, arr, table);
}
template void insert_arma_vector(std::string const &, arma::Col<double> const &,
                                 toml::table &);
template void insert_arma_vector(std::string const &,
                                 arma::Col<arma::sword> const &, toml::table &);

template <typename T>
void insert_arma_matrix(std::string const &key, arma::Mat<T> const &mat,
                        toml::table &table) {
  toml::array arr;
  for (arma::uword i = 0; i < mat.n_rows; ++i) {
    toml::array row;
    for (arma::uword j = 0; j < mat.n_cols; ++j) {
      row.push_back(mat(i, j));
    }
    arr.push_back(row);
  }
  insert(key, arr, table);
}

template <>
void insert_arma_matrix(std::string const &key,
                        arma::Mat<arma::uword> const &mat, toml::table &table) {
  using namespace arma;
  toml::array arr;
  for (arma::uword i = 0; i < mat.n_rows; ++i) {
    toml::array row;
    for (arma::uword j = 0; j < mat.n_cols; ++j) {
      row.push_back((int64_t)mat(i, j));
    }
    arr.push_back(row);
  }
  insert(key, arr, table);
}

template <>
void insert_arma_matrix(std::string const &key, arma::Mat<complex> const &mat,
                        toml::table &table) {
  using namespace arma;
  toml::array arr;
  for (arma::uword i = 0; i < mat.n_rows; ++i) {
    toml::array row;
    for (arma::uword j = 0; j < mat.n_cols; ++j) {
      row.push_back(toml::array{std::real(mat(i, j)), std::imag(mat(i, j))});
    }
    arr.push_back(row);
  }
  insert(key, arr, table);
}

template void insert_arma_matrix(std::string const &, arma::Mat<double> const &,
                                 toml::table &);
template void insert_arma_matrix(std::string const &,
                                 arma::Mat<arma::sword> const &, toml::table &);

} // namespace hydra::io
