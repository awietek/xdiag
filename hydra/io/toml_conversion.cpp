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
    return complex();
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
    return complex();
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
    return complex();
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
toml::array std_vector_to_toml_array(std::vector<T> const &value) {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(x);
  }
  return arr;
}
template <>
toml::array std_vector_to_toml_array(std::vector<uint64_t> const &value) {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back((int64_t)x);
  }
  return arr;
}

template <>
toml::array std_vector_to_toml_array(std::vector<complex> const &value) {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(toml::array{x.real(), x.imag()});
  }
  return arr;
}

template toml::array std_vector_to_toml_array(std::vector<int8_t> const &);
template toml::array std_vector_to_toml_array(std::vector<int16_t> const &);
template toml::array std_vector_to_toml_array(std::vector<int32_t> const &);
template toml::array std_vector_to_toml_array(std::vector<int64_t> const &);
template toml::array std_vector_to_toml_array(std::vector<uint8_t> const &);
template toml::array std_vector_to_toml_array(std::vector<uint16_t> const &);
template toml::array std_vector_to_toml_array(std::vector<uint32_t> const &);
template toml::array std_vector_to_toml_array(std::vector<double> const &);
template toml::array std_vector_to_toml_array(std::vector<std::string> const &);

template <typename T>
toml::array arma_vector_to_toml_array(arma::Col<T> const &value) {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(x);
  }
  return arr;
}

template <>
toml::array arma_vector_to_toml_array(arma::Col<arma::uword> const &value) {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back((int64_t)x);
  }
  return arr;
}

template <>
toml::array arma_vector_to_toml_array(arma::Col<complex> const &value) {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(toml::array{x.real(), x.imag()});
  }
  return arr;
}

template toml::array arma_vector_to_toml_array(arma::Col<double> const &);
template toml::array arma_vector_to_toml_array(arma::Col<arma::sword> const &);

template <typename T>
toml::array arma_matrix_to_toml_array(arma::Mat<T> const &mat) {
  toml::array arr;
  for (arma::uword i = 0; i < mat.n_rows; ++i) {
    toml::array row;
    for (arma::uword j = 0; j < mat.n_cols; ++j) {
      row.push_back(mat(i, j));
    }
    arr.push_back(row);
  }
  return arr;
}

template <>
toml::array arma_matrix_to_toml_array(arma::Mat<arma::uword> const &mat) {
  toml::array arr;
  for (arma::uword i = 0; i < mat.n_rows; ++i) {
    toml::array row;
    for (arma::uword j = 0; j < mat.n_cols; ++j) {
      row.push_back((int64_t)mat(i, j));
    }
    arr.push_back(row);
  }
  return arr;
}

template <>
toml::array arma_matrix_to_toml_array(arma::Mat<complex> const &mat) {
  toml::array arr;
  for (arma::uword i = 0; i < mat.n_rows; ++i) {
    toml::array row;
    for (arma::uword j = 0; j < mat.n_cols; ++j) {
      row.push_back(toml::array{std::real(mat(i, j)), std::imag(mat(i, j))});
    }
    arr.push_back(row);
  }
  return arr;
}

template toml::array arma_matrix_to_toml_array(arma::Mat<double> const &);
template toml::array arma_matrix_to_toml_array(arma::Mat<arma::sword> const &);

toml::array bond_to_toml_array(Bond const &bond) {
  toml::array array;

  // Type or Matrix
  if (bond.type_defined()) {
    array.push_back(bond.type());
  } else {
    assert(bond.matrix_defined());
    auto mat = bond.matrix();
    if (arma::norm(arma::imag(mat)) < 1e-12) {
      array.push_back(arma_matrix_to_toml_array(arma::mat(arma::real(mat))));
    } else {
      array.push_back(arma_matrix_to_toml_array(mat));
    }
  }

  // Coupling or Coupling name
  if (bond.coupling_defined()) {
    complex cpl = bond.coupling();

    if (std::abs(cpl.imag()) < 1e-12) {
      array.push_back(cpl.real());
    } else {
      array.push_back(toml::array{cpl.real(), cpl.imag()});
    }
  } else {
    assert(bond.coupling_named());
    array.push_back(bond.coupling_name());
  }

  // Sites
  if (bond.size() == 1) {
    array.push_back(bond[0]);
  } else {
    array.push_back(std_vector_to_toml_array(bond.sites()));
  }
  return array;
}

Bond toml_array_to_bond(toml::array const &array) {
  if (array.size() == 3) {
    auto type_node = array[0].value<std::string>();
    auto matrix_node = array[0].as_array();

    auto coupling_node = array[1].value<double>();
    auto coupling_node_cplx = array[1].as_array();
    auto coupling_name_node = array[1].value<std::string>();

    auto site_node = array[2].value<int>();
    auto sites_node = array[2].as_array();

    if (type_node) {
      std::string type = *type_node;
      if (coupling_node || coupling_node_cplx) {
        complex coupling = (coupling_node)
                               ? *coupling_node
                               : get_toml_value<complex>(*coupling_node_cplx);
        if (site_node) {
          int site = *site_node;
          return Bond(type, coupling, site);
        } else if (sites_node) {
          auto sites = toml_array_to_std_vector<int>(*sites_node);
          return Bond(type, coupling, sites);
        } else {
          Log.err(
              "Error parsing toml to hydra::Bond: third entry must be either a "
              "single int or vector of ints");
          return Bond();
        }

      } else if (coupling_name_node) {
        std::string coupling_name = *coupling_name_node;

        if (site_node) {
          int site = *site_node;
          return Bond(type, coupling_name, site);
        } else if (sites_node) {
          auto sites = toml_array_to_std_vector<int>(*sites_node);
          return Bond(type, coupling_name, sites);
        } else {
          Log.err(
              "Error parsing toml to hydra::Bond: third entry must be either a "
              "single int or vector of ints");
          return Bond();
        }

      } else {
        Log.err(
            "Error parsing toml to hydra::Bond: second entry must be either a "
            "string or a (complex) number");
        return Bond();
      }

    } else if (matrix_node) {
      arma::cx_mat matrix = toml_array_to_arma_matrix<complex>(*matrix_node);

      if (coupling_node || coupling_node_cplx) {
        complex coupling = (coupling_node)
                               ? *coupling_node
                               : get_toml_value<complex>(*coupling_node_cplx);

        if (site_node) {
          int site = *site_node;
          return Bond(matrix, coupling, site);
        } else if (sites_node) {
          auto sites = toml_array_to_std_vector<int>(*sites_node);
          return Bond(matrix, coupling, sites);
        } else {
          Log.err(
              "Error parsing toml to hydra::Bond: third entry must be either a "
              "single int or vector of ints");
          return Bond();
        }

      } else if (coupling_name_node) {
        std::string coupling_name = *coupling_name_node;

        if (site_node) {
          int site = *site_node;
          return Bond(matrix, coupling_name, site);
        } else if (sites_node) {
          auto sites = toml_array_to_std_vector<int>(*sites_node);
          return Bond(matrix, coupling_name, sites);
        } else {
          Log.err(
              "Error parsing toml to hydra::Bond: third entry must be either a "
              "single int or vector of ints");
          return Bond();
        }
      } else {
        Log.err(
            "Error parsing toml to hydra::Bond: second entry must be either a "
            "string or a (complex) number");
        return Bond();
      }

    } else {
      Log.err("Error parsing toml to hydra::Bond: first entry must be either a "
              "string or a matrix");
      return Bond();
    }
  } else {
    Log.err("Error parsing toml to hydra::Bond: bond is not an array of "
            "length 3");
    return Bond();
  }
}

BondList toml_array_to_bond_list(toml::array const &array) {
  BondList bonds;
  for (std::size_t i = 0; i < array.size(); ++i) {
    auto bond_array = array[i].as_array();
    if (bond_array) {
      bonds << toml_array_to_bond(*bond_array);
    } else {
      Log.err("Error parsing toml to hydra::BondList: entry {} is not a "
              "toml::array",
              i);
    }
  }
  return bonds;
}

BondList toml_table_to_bond_list(toml::table const &table) {
  BondList bonds;
  auto interactions_opt = table["Interactions"].as_array();
  if (interactions_opt) {
    bonds = toml_array_to_bond_list(*interactions_opt);
  }

  auto couplings_opt = table["Couplings"].as_table();
  if (couplings_opt) {
    toml::table couplings = *couplings_opt;
    for (auto [key_node, value_node] : couplings) {
      std::string key(key_node.str());
      auto value = get_toml_value<complex>(value_node);
      bonds.set_coupling(key, value);
    }
  }

  auto matrices_opt = table["Matrices"].as_table();
  if (matrices_opt) {
    toml::table matrices = *matrices_opt;
    for (auto [key_node, value_node] : matrices) {
      std::string key(key_node.str());
      auto matrix_arr_opt = value_node.as_array();
      if (matrix_arr_opt) {
        auto matrix = toml_array_to_arma_matrix<complex>(*matrix_arr_opt);
        bonds.set_matrix(key, matrix);
      } else {
        Log.err("Error parsing matrix of BondList: matrix is not a toml array");
      }
    }
  }
  return bonds;
}

toml::array bond_list_to_toml_array(BondList const &bonds) {
  toml::array array;
  for (auto &&bond : bonds) {
    array.push_back(bond_to_toml_array(bond));
  }
  return array;
}

toml::table bond_list_to_toml_table(BondList const &bonds) {
  toml::table table;
  table.insert_or_assign("Interactions", bond_list_to_toml_array(bonds));

  if (bonds.couplings().size() > 0) {
    toml::table couplings;
    for (auto &&[name, coupling] : bonds.couplings()) {
      if (std::abs(coupling.imag()) < 1e-12) {
        couplings.insert_or_assign(name, coupling.real());
      } else {
        couplings.insert_or_assign(
            name, toml::array{coupling.real(), coupling.imag()});
      }
    }
    table.insert_or_assign("Couplings", couplings);
  }

  if (bonds.matrices().size() > 0) {
    toml::table matrices;
    for (auto &&[name, matrix] : bonds.matrices()) {
      if (arma::norm(arma::imag(matrix)) < 1e-12) {
        matrices.insert_or_assign(
            name, arma_matrix_to_toml_array(arma::mat(arma::real(matrix))));
      } else {
        matrices.insert_or_assign(name, arma_matrix_to_toml_array(matrix));
      }
    }
    table.insert_or_assign("Matrices", matrices);
  }

  return table;
}

} // namespace hydra::io
