#include "toml_conversion.hpp"

#include <cassert>

#include <xdiag/common.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::io {

template <typename T> T get_toml_value(toml_const_node_view const &node) try {
  auto val = node.value<T>();
  if (val) {
    return T(*val);
  } else {
    XDIAG_THROW("Error parsing toml: node does not contain a value!");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
complex get_toml_value<complex>(toml_const_node_view const &node) try {
  auto val = node.value<double>();
  auto arr = node.as_array();
  if (val) {
    return complex(*val);
  } else if (arr) {
    auto real_imag = toml_array_to_std_array<double, 2>(*arr);
    return complex{real_imag[0], real_imag[1]};
  } else {
    XDIAG_THROW("Error parsing toml: failed trying to parse complex number!");
    return complex();
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
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

template <typename T> T get_toml_value(toml_node_view const &node) try {
  auto val = node.value<T>();
  if (val) {
    return T(*val);
  } else {
    XDIAG_THROW("Error parsing toml: node does not contain a value!");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> complex get_toml_value<complex>(toml_node_view const &node) try {
  auto val = node.value<double>();
  auto arr = node.as_array();
  if (val) {
    return complex(*val);
  } else if (arr) {
    auto real_imag = toml_array_to_std_array<double, 2>(*arr);
    return complex{real_imag[0], real_imag[1]};
  } else {
    XDIAG_THROW("Error parsing toml: failed trying to parse complex number!");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
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

template <typename T> T get_toml_value(toml::node const &node) try {
  auto val = node.value<T>();
  if (val) {
    return T(*val);
  } else {
    XDIAG_THROW("Error parsing toml: node does not contain a value!");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> complex get_toml_value<complex>(toml::node const &node) try {
  auto val = node.value<double>();
  auto arr = node.as_array();
  if (val) {
    return complex(*val);
  } else if (arr) {
    auto real_imag = toml_array_to_std_array<double, 2>(*arr);
    return complex{real_imag[0], real_imag[1]};
  } else {
    XDIAG_THROW("Error parsing toml: failed trying to parse complex number!");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
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

toml::array get_toml_array(toml_const_node_view const &node) try {
  auto arr = node.as_array();
  if (arr) {
    return *arr;
  } else {
    XDIAG_THROW("Error parsing toml: node does not contain an array!");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

toml::array get_toml_array(toml_node_view const &node) try {
  auto arr = node.as_array();
  if (arr) {
    return *arr;
  } else {
    XDIAG_THROW("Error parsing toml: node does not contain an array!");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

toml::array get_toml_array(toml::node const &node) try {
  auto arr = node.as_array();
  if (arr) {
    return *arr;
  } else {
    XDIAG_THROW("Error parsing toml: node does not contain an array!");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename T, std::size_t size>
std::array<T, size> toml_array_to_std_array(toml::array const &toml_array) try {
  if (size != toml_array.size()) {
    XDIAG_THROW(fmt::format(
        "Error reading TOML array to std::array: toml::array size ({}) "
        "exceeds std::array size ({}) !",
        toml_array.size(), size));
  }

  std::array<T, size> array;
  for (std::size_t i = 0; i < size; ++i) {
    array[i] = get_toml_value<T>(toml_array[i]);
  }
  return array;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template std::array<double, 2>
toml_array_to_std_array<double, 2>(toml::array const &);

template <typename T>
std::vector<T> toml_array_to_std_vector(toml::array const &toml_array) try {
  std::size_t size = toml_array.size();
  std::vector<T> vector(size);
  for (std::size_t i = 0; i < size; ++i) {
    vector[i] = get_toml_value<T>(toml_array[i]);
  }
  return vector;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
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
arma::Col<T> toml_array_to_arma_vector(toml::array const &toml_array) try {
  std::size_t size = toml_array.size();
  arma::Col<T> vector(size);
  for (arma::uword i = 0; i < size; ++i) {
    vector[i] = get_toml_value<T>(toml_array[i]);
  }
  return vector;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
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
arma::Mat<T> toml_array_to_arma_matrix(toml::array const &toml_array) try {
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
          XDIAG_THROW("Error reading TOML array to arma::Matrix: not all rows "
                      "have same length");
        }

        for (std::size_t j = 0; j < n; ++j) {
          vector[i + m * j] = get_toml_value<T>(row_entries[j]);
        }

      } else {
        XDIAG_THROW(fmt::format(
            "Error reading TOML array to arma::Matrix: Cannot parse row {}!",
            i));
      }
    } // for (int i = 0; i < m; ++i)
    return arma::Mat<T>(vector.data(), m, n);
  } else {
    return arma::Mat<T>();
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
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
toml::array std_vector_to_toml_array(std::vector<T> const &value) try {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(x);
  }
  return arr;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
toml::array std_vector_to_toml_array(std::vector<uint64_t> const &value) try {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back((int64_t)x);
  }
  return arr;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
toml::array std_vector_to_toml_array(std::vector<complex> const &value) try {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(toml::array{x.real(), x.imag()});
  }
  return arr;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
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
toml::array arma_vector_to_toml_array(arma::Col<T> const &value) try {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(x);
  }
  return arr;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
toml::array arma_vector_to_toml_array(arma::Col<arma::uword> const &value) try {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back((int64_t)x);
  }
  return arr;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
toml::array arma_vector_to_toml_array(arma::Col<complex> const &value) try {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(toml::array{x.real(), x.imag()});
  }
  return arr;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template toml::array arma_vector_to_toml_array(arma::Col<double> const &);
template toml::array arma_vector_to_toml_array(arma::Col<arma::sword> const &);

template <typename T>
toml::array arma_matrix_to_toml_array(arma::Mat<T> const &mat) try {
  toml::array arr;
  for (arma::uword i = 0; i < mat.n_rows; ++i) {
    toml::array row;
    for (arma::uword j = 0; j < mat.n_cols; ++j) {
      row.push_back(mat(i, j));
    }
    arr.push_back(row);
  }
  return arr;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
toml::array arma_matrix_to_toml_array(arma::Mat<arma::uword> const &mat) try {
  toml::array arr;
  for (arma::uword i = 0; i < mat.n_rows; ++i) {
    toml::array row;
    for (arma::uword j = 0; j < mat.n_cols; ++j) {
      row.push_back((int64_t)mat(i, j));
    }
    arr.push_back(row);
  }
  return arr;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
toml::array arma_matrix_to_toml_array(arma::Mat<complex> const &mat) try {
  toml::array arr;
  for (arma::uword i = 0; i < mat.n_rows; ++i) {
    toml::array row;
    for (arma::uword j = 0; j < mat.n_cols; ++j) {
      row.push_back(toml::array{std::real(mat(i, j)), std::imag(mat(i, j))});
    }
    arr.push_back(row);
  }
  return arr;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return toml::array();
}

template toml::array arma_matrix_to_toml_array(arma::Mat<double> const &);
template toml::array arma_matrix_to_toml_array(arma::Mat<arma::sword> const &);

toml::array bond_to_toml_array(Bond const &bond) try {
  toml::array array;

  // Type
  array.push_back(bond.type());

  // Coupling
  if (bond.coupling().is<std::string>()) {
    array.push_back(bond.coupling().as<std::string>());
  } else if (bond.coupling().is<double>()) {
    array.push_back(bond.coupling().as<double>());
  } else if (bond.coupling().is<complex>()) {
    complex cpl = bond.coupling().as<complex>();
    array.push_back(toml::array{cpl.real(), cpl.imag()});
  } else if (bond.coupling().is<arma::mat>()) {
    arma::mat mat = bond.coupling().as<arma::mat>();
    array.push_back(arma_matrix_to_toml_array(mat));
  } else if (bond.coupling().is<arma::cx_mat>()) {
    arma::cx_mat mat = bond.coupling().as<arma::cx_mat>();
    array.push_back(arma_matrix_to_toml_array(mat));
  } else {
    XDIAG_THROW("Unable to convert Bond to toml array ");
  }

  // Sites
  for (int i = 0; i < bond.size(); ++i) {
    array.push_back(bond[i]);
  }

  return array;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Bond toml_array_to_bond(toml::array const &array) try {
  if (array.size() < 3) {
    XDIAG_THROW(
        "Error parsing toml to xdiag::Bond: toml array must have exactly three "
        "entries.\n"
        "1) A string defining the type of the bond.\n"
        "2) A coupling that can either be a string, a "
        "real/complex number, or a real/complex matrix.\n"
        "3)  thesites of the bond must be either a single integer (single site "
        "bond) or a 1D array of integers.")
  }

  // First get the type
  std::string type;
  auto type_node = array[0].value<std::string>();
  if (type_node) {
    type = *type_node;
  } else {
    XDIAG_THROW(
        "Error parsing toml to xdiag::Bond: first entry must be a string ");
  }

  // then get the sites
  std::vector<int64_t> sites;
  for (int i = 2; i < array.size(); ++i) {
    auto site_node = array[i].value<int64_t>();
    if (site_node) {
      sites.push_back(*site_node);

    } else {
      XDIAG_THROW("Error parsing toml to xdiag::Bond: third and onward entries "
                  "must be integers");
    }
  }
  // then get the coupling
  auto coupling_node_string = array[1].value<std::string>();
  auto coupling_node_double = array[1].value<double>();
  auto coupling_node_array = array[1].as_array();
  if (coupling_node_string) {
    std::string coupling = *coupling_node_string;
    return Bond(type, coupling, sites);
  } else if (coupling_node_double) {
    double coupling = *coupling_node_double;
    return Bond(type, coupling, sites);
  } else if (coupling_node_array) {
    auto coupling_array = *coupling_node_array;
    if (coupling_array.size() == 0) {
      XDIAG_THROW("Error parsing toml to xdiag::Bond: an array was handed to "
                  "the second entry defining the coupling. However, it is "
                  "found to be empty");
    }

    auto node_real = coupling_array[1].value<double>();
    auto node_array = array[1].as_array();
    // coupling should be a complex number
    if (node_real) {
      complex coupling = get_toml_value<complex>(coupling_array);
      return Bond(type, coupling, sites);
    }
    // Coupling is either a real or complex matrix
    else if (node_array) {
      auto mat_array = *node_array;
      if (mat_array.size() == 0) {
        XDIAG_THROW("Error parsing toml to xdiag::Bond: an array was handed to "
                    "the second entry defining the coupling. However, it is "
                    "found to be empty");
      }
      auto mat2_array = mat_array[0].as_array();
      if (mat2_array) {
        if (mat2_array->size() == 0) {
          XDIAG_THROW(
              "Error parsing toml to xdiag::Bond: an array was handed to "
              "the second entry defining the coupling. However, it is "
              "found to be empty");
        }
        auto real_array = (*mat2_array)[0].value<double>();
        auto cplx_array = (*mat2_array)[0].as_array();

        if (real_array) {
          arma::mat coupling = toml_array_to_arma_matrix<double>(mat_array);
          return Bond(type, coupling, sites);
        } else if (cplx_array) {
          arma::cx_mat coupling = toml_array_to_arma_matrix<complex>(mat_array);
          return Bond(type, coupling, sites);
        } else {
          XDIAG_THROW("Error parsing toml to xdiag::Bond: Invalid form of the "
                      "coupling.");
        }
      } else {
        XDIAG_THROW("Error parsing toml to xdiag::Bond: Invalid form of the "
                    "coupling.");
      }
    } else {
      XDIAG_THROW("Error parsing toml to xdiag::Bond: Invalid form of the "
                  "coupling.");
    }
  } else {
    XDIAG_THROW("Error parsing toml to xdiag::Bond: Invalid form of the "
                "coupling.");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

BondList toml_array_to_bond_list(toml::array const &array) try {
  BondList bonds;
  for (std::size_t i = 0; i < array.size(); ++i) {
    auto bond_array = array[i].as_array();
    if (bond_array) {
      bonds += toml_array_to_bond(*bond_array);
    } else {
      XDIAG_THROW(fmt::format(
          "Error parsing toml to xdiag::BondList: entry {} is not a "
          "toml::array",
          i));
    }
  }
  return bonds;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

BondList toml_table_to_bond_list(toml::table const &table) try {
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

      auto val = value_node.value<double>();
      if (val) {
        auto value = get_toml_value<double>(value_node);
        bonds[key] = value;
      } else {
        if (value_node.as_array()->size() == 0) {
          XDIAG_THROW("Invalid format of coupling in BondList")
        }

        auto val = *value_node.as_array();
        auto val_real = val[0].value<double>();
        auto val_array = val[0].as_array();
        if (val_real) {
          auto value = get_toml_value<complex>(value_node);
          bonds[key] = value;
        }
        if (val_array) {
          auto val = *(val_array->as_array());

          auto val2_real = val[0].value<double>();
          auto val2_arr = val[0].as_array();
          if (val2_real) {
            auto value =
                toml_array_to_arma_matrix<double>(*value_node.as_array());
            bonds[key] = value;
          } else if (val2_arr) {
            auto value =
                toml_array_to_arma_matrix<complex>(*value_node.as_array());
            bonds[key] = value;

          } else {
            XDIAG_THROW("Invalid format of coupling in BondList")
          }
        }
      }
    }
  }

  return bonds;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

toml::array bond_list_to_toml_array(BondList const &bonds) try {
  toml::array array;
  for (auto &&bond : bonds) {
    array.push_back(bond_to_toml_array(bond));
  }
  return array;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

toml::table bond_list_to_toml_table(BondList const &bonds) try {
  toml::table table;
  table.insert_or_assign("Interactions", bond_list_to_toml_array(bonds));

  if (bonds.couplings().size() > 0) {
    toml::table couplings;
    for (auto name : bonds.couplings()) {
      Coupling cpl = bonds[name];
      if (cpl.is<std::string>()) {
        couplings.insert_or_assign(name, cpl.as<std::string>());
      } else if (cpl.is<double>()) {
        couplings.insert_or_assign(name, cpl.as<double>());
      } else if (cpl.is<complex>()) {
        complex c = cpl.as<complex>();
        couplings.insert_or_assign(name, toml::array{c.real(), c.imag()});
      } else if (cpl.is<arma::mat>()) {
        couplings.insert_or_assign(
            name, arma_matrix_to_toml_array(cpl.as<arma::mat>()));
      } else if (cpl.is<arma::cx_mat>()) {
        couplings.insert_or_assign(
            name, arma_matrix_to_toml_array(cpl.as<arma::cx_mat>()));
      } else {
        XDIAG_THROW("Unable to convert BondList to toml");
      }
      table.insert_or_assign("Couplings", couplings);
    }
  }
  return table;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::io
