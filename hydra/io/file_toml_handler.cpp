#include "file_toml_handler.h"

#include <extern/armadillo/armadillo>
#include <extern/toml++/toml.h>

#include <array>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <vector>

#include <hydra/common.h>
#include <hydra/utils/logger.h>

namespace hydra::io {

FileTomlHandler::FileTomlHandler(std::string key, toml::table &table)
    : key_(key), table_(table) {}

template <typename T>
T retrieve_plain(std::string const &key, toml::table const &table) {
  auto val = table.at_path(key).value<T>();
  if (val) {
    return T(*val);
  } else {
    Log.err("Error parsing toml: key \"{}\" not found!", key);
    return T();
  }
}

template <class T, int size>
std::array<T, size> retrieve_std_array(std::string const &key,
                                       toml::table const &table) {
  std::array<T, size> array;

  auto toml_array_opt = table.at_path(key).as_array();
  if (toml_array_opt) {
    toml::array toml_array = *toml_array_opt;
    if (size != toml_array.size()) {
      Log.err("Error reading TOML array to std::array: toml::array size ({}) "
              "exceeds std::array size ({}) !",
              toml_array.size(), size);
    }

    for (int i = 0; i < size; ++i) {
      auto val_opt = toml_array[i].value<T>();
      if (val_opt) {
        array[i] = *val_opt;
      } else {
        Log.err("Error reading TOML array to std::array: cannot properly "
                "convert the type of the entry!");
      }
    }
  } else {
    Log.err("Error reading TOML array to std::array: key \"{}\" not found!",
            key);
  }
  return array;
}

template <typename T>
std::vector<T> retrieve_std_vector(std::string const &key,
                                   toml::table const &table) {
  std::vector<T> vector;

  auto toml_array_opt = table.at_path(key).as_array();
  if (toml_array_opt) {
    toml::array toml_array = *toml_array_opt;
    int size = toml_array.size();

    std::vector<T> vector(size);

    for (int i = 0; i < size; ++i) {
      auto val_opt = toml_array[i].value<T>();
      if (val_opt) {
        vector[i] = *val_opt;
      } else {
        Log.err("Error reading TOML array to std::vector: cannot properly "
                "convert the type of the entry!");
      }
    }

    return vector;
  } else {
    Log.err("Error reading TOML array to std::vector: key \"{}\" not found!",
            key);
    return std::vector<T>();
  }
}

template <>
std::vector<complex> retrieve_std_vector<complex>(std::string const &key,
                                                  toml::table const &table) {
  std::vector<complex> vector;
  auto toml_array_opt = table.at_path(key).as_array();
  if (toml_array_opt) {
    toml::array toml_array = *toml_array_opt;
    int size = toml_array.size();

    std::vector<complex> vector(size);
    for (int i = 0; i < size; ++i) {
      auto cplx_opt = toml_array[i].as_array();
      if (cplx_opt) {
        toml::array cplx_array = *cplx_opt;
        if (cplx_array.size() == 2) {
          std::array<double, 2> cplx_std_array = {0., 0.};
          for (int j = 0; j < 2; ++j) {
            auto val_opt = cplx_array[j].value<double>();
            if (val_opt) {
              cplx_std_array[j] = *val_opt;
            } else {
              Log.err("Error reading TOML array to std::vector<complex>: "
                      "cannot properly convert an entry to double!");
            }
          }
          vector[i] = complex{cplx_std_array[0], cplx_std_array[1]};
        } else {
          Log.err("Error reading TOML array to std::vector<complex>: complex "
                  "numbers need to be specified as length 2 toml::arrays!");
        }

      } else {
        Log.err("Error reading TOML array to std::vector<complex>: entries of "
                "the main toml::array are not toml::arrays");
      }
    }

    return vector;
  } else {
    Log.err("Error reading TOML array to std::vector: key \"{}\" not found!",
            key);
    return std::vector<complex>();
  }
}

// Integer types
template <> bool FileTomlHandler::as<bool>() {
  return retrieve_plain<bool>(key_, table_);
}

template <> int8_t FileTomlHandler::as<int8_t>() {
  return retrieve_plain<int8_t>(key_, table_);
}

template <> int16_t FileTomlHandler::as<int16_t>() {
  return retrieve_plain<int16_t>(key_, table_);
}

template <> int32_t FileTomlHandler::as<int32_t>() {
  return retrieve_plain<int32_t>(key_, table_);
}

template <> int64_t FileTomlHandler::as<int64_t>() {
  return retrieve_plain<int64_t>(key_, table_);
}

template <> uint8_t FileTomlHandler::as<uint8_t>() {
  return retrieve_plain<uint8_t>(key_, table_);
}

template <> uint16_t FileTomlHandler::as<uint16_t>() {
  return retrieve_plain<uint16_t>(key_, table_);
}

template <> uint32_t FileTomlHandler::as<uint32_t>() {
  return retrieve_plain<uint32_t>(key_, table_);
}

template <> uint64_t FileTomlHandler::as<uint64_t>() {
  return retrieve_plain<uint64_t>(key_, table_);
}

// Floating points
template <> double FileTomlHandler::as<double>() {
  return retrieve_plain<double>(key_, table_);
}

template <> std::complex<double> FileTomlHandler::as<std::complex<double>>() {
  auto arr = retrieve_std_array<double, 2>(key_, table_);
  return std::complex<double>(arr[0], arr[1]);
}

// strings
template <> std::string FileTomlHandler::as<std::string>() {
  return retrieve_plain<std::string>(key_, table_);
}

// std::vectors
template <> std::vector<int8_t> FileTomlHandler::as<std::vector<int8_t>>() {
  return retrieve_std_vector<int8_t>(key_, table_);
}
template <> std::vector<int16_t> FileTomlHandler::as<std::vector<int16_t>>() {
  return retrieve_std_vector<int16_t>(key_, table_);
}
template <> std::vector<int32_t> FileTomlHandler::as<std::vector<int32_t>>() {
  return retrieve_std_vector<int32_t>(key_, table_);
}
template <> std::vector<int64_t> FileTomlHandler::as<std::vector<int64_t>>() {
  return retrieve_std_vector<int64_t>(key_, table_);
}
template <> std::vector<uint8_t> FileTomlHandler::as<std::vector<uint8_t>>() {
  return retrieve_std_vector<uint8_t>(key_, table_);
}
template <> std::vector<uint16_t> FileTomlHandler::as<std::vector<uint16_t>>() {
  return retrieve_std_vector<uint16_t>(key_, table_);
}
template <> std::vector<uint32_t> FileTomlHandler::as<std::vector<uint32_t>>() {
  return retrieve_std_vector<uint32_t>(key_, table_);
}
template <> std::vector<uint64_t> FileTomlHandler::as<std::vector<uint64_t>>() {
  return retrieve_std_vector<uint64_t>(key_, table_);
}

template <> std::vector<double> FileTomlHandler::as<std::vector<double>>() {
  return retrieve_std_vector<double>(key_, table_);
}
template <> std::vector<complex> FileTomlHandler::as<std::vector<complex>>() {
  return retrieve_std_vector<complex>(key_, table_);
}
template <>
std::vector<std::string> FileTomlHandler::as<std::vector<std::string>>() {
  return retrieve_std_vector<std::string>(key_, table_);
}

// Armadillo stuff
template <> arma::vec FileTomlHandler::as<arma::vec>() {
  return arma::vec(retrieve_std_vector<double>(key_, table_));
}

template <> arma::cx_vec FileTomlHandler::as<arma::cx_vec>() {
  return arma::cx_vec(retrieve_std_vector<complex>(key_, table_));
}

template <> arma::ivec FileTomlHandler::as<arma::ivec>() {
  return arma::ivec(retrieve_std_vector<arma::sword>(key_, table_));
}

template <> arma::uvec FileTomlHandler::as<arma::uvec>() {
  return arma::uvec(retrieve_std_vector<arma::uword>(key_, table_));
}

template <typename T>
arma::Mat<T> retrieve_arma_matrix(std::string const &key,
                                  toml::table const &table) {
  auto rows_opt = table.at_path(key).as_array();
  if (rows_opt) {
    toml::array rows = *rows_opt;
    int m = rows.size();
    int n = 0;
    if (m > 0) {
      std::vector<T> vector;
      for (int i = 0; i < m; ++i) {
        auto row_entries_opt = rows[i].as_array();
        if (row_entries_opt) {
          toml::array row_entries = *row_entries_opt;
          if (i == 0) {
            n = row_entries.size();
            vector.resize(m * n);
          }
          if ((int)row_entries.size() != n) {
            Log.err("Error reading TOML array to arma::Matrix: not all rows "
                    "have same length");
          }
          for (int j = 0; j < n; ++j) {
            auto val_opt = row_entries[j].value<T>();
            if (val_opt) {
              vector[i + m * j] = *val_opt;
            } else {
              Log.err(
                  "Error reading TOML array to arma::Matrix: cannot properly "
                  "convert the type of the entry!");
            }
          }
        } else {
          Log.err(
              "Error reading TOML array to arma::Matrix: Cannot parse row {}!",
              i);
        }
      }
      return arma::Mat<T>(vector.data(), m, n);
    } else {
      return arma::Mat<T>();
    }
  } else {
    Log.err("Error parsing toml: key \"{}\" not found!", key);
    return arma::Mat<T>();
  }
}

template <>
arma::Mat<complex> retrieve_arma_matrix<complex>(std::string const &key,
                                                 toml::table const &table) {
  auto rows_opt = table.at_path(key).as_array();
  if (rows_opt) {
    toml::array rows = *rows_opt;
    int m = rows.size();
    int n = 0;
    if (m > 0) {
      std::vector<complex> vector;
      for (int i = 0; i < m; ++i) {
        auto row_entries_opt = rows[i].as_array();
        if (row_entries_opt) {
          toml::array row_entries = *row_entries_opt;
          if (i == 0) {
            n = row_entries.size();
            vector.resize(m * n);
          }
          if ((int)row_entries.size() != n) {
            Log.err("Error reading TOML array to arma::Matrix: not all rows "
                    "have same length");
          }
          for (int j = 0; j < n; ++j) {
            auto val_opt = row_entries[j].value<double>();
            auto arr_opt = row_entries[j].as_array();
            if (val_opt) {
              vector[i + m * j] = complex(*val_opt);
            } else if (arr_opt) {
              toml::array cplx_array = *arr_opt;
              if (cplx_array.size() == 2) {
                std::array<double, 2> cplx_std_array = {0., 0.};
                for (int k = 0; k < 2; ++k) {
                  auto c_opt = cplx_array[k].value<double>();
                  if (c_opt) {
                    cplx_std_array[k] = *c_opt;
                  } else {
                    Log.err(
                        "Error reading TOML array to arma::Matrix<complex>: "
                        "cannot properly convert an entry to double!");
                  }
                }
                vector[i + m * j] =
                    complex{cplx_std_array[0], cplx_std_array[1]};
              } else {
                Log.err(
                    "Error reading TOML array to arma::Matrix<complex>: "
                    "complex "
                    "numbers need to be specified as length 2 toml::arrays!");
              }
            } else {
              Log.err(
                  "Error reading TOML array to arma::Matrix: cannot properly "
                  "convert the type of the entry!");
            }
          }

        } else {
          Log.err("Error reading TOML array to arma::Matrix: Cannot parse "
                  "row {}!",
                  i);
        }
      }
      return arma::Mat<complex>(vector.data(), m, n);
    } else {
      return arma::Mat<complex>();
    }
  } else {
    Log.err("Error parsing toml: key \"{}\" not found!", key);
    return arma::Mat<complex>();
  }
}

template <> arma::mat FileTomlHandler::as<arma::mat>() {
  return arma::mat(retrieve_arma_matrix<double>(key_, table_));
}

template <> arma::cx_mat FileTomlHandler::as<arma::cx_mat>() {
  return arma::cx_mat(retrieve_arma_matrix<complex>(key_, table_));
}

template <> arma::imat FileTomlHandler::as<arma::imat>() {
  return arma::imat(retrieve_arma_matrix<arma::sword>(key_, table_));
}

template <> arma::umat FileTomlHandler::as<arma::umat>() {
  return arma::umat(retrieve_arma_matrix<arma::uword>(key_, table_));
}

//////////////////////////////////////////////////////////////////
// operator=

template <typename T> void FileTomlHandler::operator=(T const &value) {
  table_.insert_or_assign(key_, value);
}
template void FileTomlHandler::operator=<std::string>(std::string const &value);
template void FileTomlHandler::operator=<int8_t>(int8_t const &value);
template void FileTomlHandler::operator=<int16_t>(int16_t const &value);
template void FileTomlHandler::operator=<int32_t>(int32_t const &value);
template void FileTomlHandler::operator=<int64_t>(int64_t const &value);
template void FileTomlHandler::operator=<uint8_t>(uint8_t const &value);
template void FileTomlHandler::operator=<uint16_t>(uint16_t const &value);
template void FileTomlHandler::operator=<uint32_t>(uint32_t const &value);
template void FileTomlHandler::operator=<double>(double const &value);

template <> void FileTomlHandler::operator=<uint64_t>(uint64_t const &value) {
  table_.insert_or_assign(key_, (int64_t)value);
}

template <> void FileTomlHandler::operator=(complex const &value) {
  table_.insert_or_assign(key_,
                          toml::array{std::real(value), std::imag(value)});
}

template <typename T>
void insert_vector(std::string const &key, toml::table &table,
                   std::vector<T> const &vec) {
  toml::array arr;
  for (auto x : vec) {
    arr.push_back(x);
  }
  table.insert_or_assign(key, arr);
}

template <>
void insert_vector(std::string const &key, toml::table &table,
                   std::vector<uint64_t> const &vec) {
  toml::array arr;
  for (auto x : vec) {
    arr.push_back((int64_t)x);
  }
  table.insert_or_assign(key, arr);
}

template <>
void insert_vector(std::string const &key, toml::table &table,
                   std::vector<complex> const &vec) {
  toml::array arr;
  for (auto x : vec) {
    arr.push_back(toml::array{std::real(x), std::imag(x)});
  }
  table.insert_or_assign(key, arr);
}

template <>
void FileTomlHandler::operator=
    <std::vector<int8_t>>(std::vector<int8_t> const &value) {
  insert_vector(key_, table_, value);
}
template <>
void FileTomlHandler::operator=
    <std::vector<int16_t>>(std::vector<int16_t> const &value) {
  insert_vector(key_, table_, value);
}
template <>
void FileTomlHandler::operator=
    <std::vector<int32_t>>(std::vector<int32_t> const &value) {
  insert_vector(key_, table_, value);
}
template <>
void FileTomlHandler::operator=
    <std::vector<int64_t>>(std::vector<int64_t> const &value) {
  insert_vector(key_, table_, value);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint8_t>>(std::vector<uint8_t> const &value) {
  insert_vector(key_, table_, value);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint16_t>>(std::vector<uint16_t> const &value) {
  insert_vector(key_, table_, value);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint32_t>>(std::vector<uint32_t> const &value) {
  insert_vector(key_, table_, value);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint64_t>>(std::vector<uint64_t> const &value) {
  insert_vector(key_, table_, value);
}
template <>
void FileTomlHandler::operator=
    <std::vector<double>>(std::vector<double> const &value) {
  insert_vector(key_, table_, value);
}
template <>
void FileTomlHandler::operator=
    <std::vector<complex>>(std::vector<complex> const &value) {
  insert_vector(key_, table_, value);
}

template <>
void FileTomlHandler::operator=
    <std::vector<std::string>>(std::vector<std::string> const &value) {
  insert_vector(key_, table_, value);
}

template <> void FileTomlHandler::operator=<arma::vec>(arma::vec const &value) {
  insert_vector(key_, table_, arma::conv_to<std::vector<double>>::from(value));
}
template <>
void FileTomlHandler::operator=<arma::cx_vec>(arma::cx_vec const &value) {
  insert_vector(key_, table_, arma::conv_to<std::vector<complex>>::from(value));
}
template <>
void FileTomlHandler::operator=<arma::ivec>(arma::ivec const &value) {
  insert_vector(key_, table_, arma::conv_to<std::vector<int64_t>>::from(value));
}
template <>
void FileTomlHandler::operator=<arma::uvec>(arma::uvec const &value) {
  insert_vector(key_, table_,
                arma::conv_to<std::vector<uint64_t>>::from(value));
}

template <typename T>
void insert_matrix(std::string const &key, toml::table &table,
                   arma::Mat<T> const &mat) {
  using namespace arma;
  toml::array arr;
  for (uword i = 0; i < mat.n_rows; ++i) {
    toml::array row;
    for (uword j = 0; j < mat.n_cols; ++j) {
      row.push_back(mat(i, j));
    }
    arr.push_back(row);
  }
  table.insert_or_assign(key, arr);
}

template <>
void insert_matrix(std::string const &key, toml::table &table,
                   arma::Mat<arma::uword> const &mat) {
  using namespace arma;
  toml::array arr;
  for (uword i = 0; i < mat.n_rows; ++i) {
    toml::array row;
    for (uword j = 0; j < mat.n_cols; ++j) {
      row.push_back((int64_t)mat(i, j));
    }
    arr.push_back(row);
  }
  table.insert_or_assign(key, arr);
}

template <>
void insert_matrix(std::string const &key, toml::table &table,
                   arma::Mat<complex> const &mat) {
  using namespace arma;
  toml::array arr;
  for (uword i = 0; i < mat.n_rows; ++i) {
    toml::array row;
    for (uword j = 0; j < mat.n_cols; ++j) {
      row.push_back(toml::array{std::real(mat(i, j)), std::imag(mat(i, j))});
    }
    arr.push_back(row);
  }
  table.insert_or_assign(key, arr);
}

template <> void FileTomlHandler::operator=<arma::mat>(arma::mat const &value) {
  insert_matrix(key_, table_, value);
}
template <>
void FileTomlHandler::operator=<arma::cx_mat>(arma::cx_mat const &value) {
  insert_matrix(key_, table_, value);
}
template <>
void FileTomlHandler::operator=<arma::imat>(arma::imat const &value) {
  insert_matrix(key_, table_, value);
}
template <>
void FileTomlHandler::operator=<arma::umat>(arma::umat const &value) {
  insert_matrix(key_, table_, value);
}

} // namespace hydra::io
