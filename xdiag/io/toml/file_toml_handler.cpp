#include "file_toml_handler.hpp"

#include <array>
#include <complex>
#include <cstdint>
#include <vector>

#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/extern/toml++/toml.hpp>

#include <xdiag/io/toml/arma_matrix.hpp>
#include <xdiag/io/toml/arma_vector.hpp>
#include <xdiag/io/toml/operators.hpp>
#include <xdiag/io/toml/std_vector.hpp>
#include <xdiag/io/toml/toml_conversion.hpp>
#include <xdiag/io/toml/value.hpp>

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/type_string.hpp>

namespace xdiag::io {

FileTomlHandler::FileTomlHandler(std::string key, toml::table &table)
    : key_(key), table_(table) {}

template <typename T>
static T as_plain(std::string key, toml::table const &table) try {
  auto node = table.at_path(key).node();
  if (node) {
    return value<T>(*node);
  } else {
    XDIAG_THROW(fmt::format("Key \"{}\" is not contained in TOML table", key));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

// Plain values
template <>
bool FileTomlHandler::as<bool>() const
    XDIAG_TRY_CATCH(return as_plain<bool>(key_, table_));

template <>
int8_t FileTomlHandler::as<int8_t>() const
    XDIAG_TRY_CATCH(return as_plain<int8_t>(key_, table_));

template <>
int16_t FileTomlHandler::as<int16_t>() const
    XDIAG_TRY_CATCH(return as_plain<int16_t>(key_, table_));

template <>
int32_t FileTomlHandler::as<int32_t>() const
    XDIAG_TRY_CATCH(return as_plain<int32_t>(key_, table_));

template <>
int64_t FileTomlHandler::as<int64_t>() const
    XDIAG_TRY_CATCH(return as_plain<int64_t>(key_, table_));

template <>
uint8_t FileTomlHandler::as<uint8_t>() const
    XDIAG_TRY_CATCH(return as_plain<uint8_t>(key_, table_));

template <>
uint16_t FileTomlHandler::as<uint16_t>() const
    XDIAG_TRY_CATCH(return as_plain<uint16_t>(key_, table_));

template <>
uint32_t FileTomlHandler::as<uint32_t>() const
    XDIAG_TRY_CATCH(return as_plain<uint32_t>(key_, table_));

template <>
uint64_t FileTomlHandler::as<uint64_t>() const
    XDIAG_TRY_CATCH(return as_plain<uint64_t>(key_, table_));

template <>
double FileTomlHandler::as<double>() const
    XDIAG_TRY_CATCH(return as_plain<double>(key_, table_));

template <>
complex FileTomlHandler::as<complex>() const
    XDIAG_TRY_CATCH(return as_plain<complex>(key_, table_));

template <>
std::string FileTomlHandler::as<std::string>() const
    XDIAG_TRY_CATCH(return as_plain<std::string>(key_, table_));

// std::vectors
template <typename T>
static std::vector<T> as_std_vector(std::string key,
                                    toml::table const &table) try {
  auto node = table.at_path(key).node();
  if (node) {
    return std_vector<T>(*node);
  } else {
    XDIAG_THROW(fmt::format("Key \"{}\" is not contained in TOML table", key));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
std::vector<int8_t> FileTomlHandler::as<std::vector<int8_t>>() const
    XDIAG_TRY_CATCH(return as_std_vector<int8_t>(key_, table_));

template <>
std::vector<int16_t> FileTomlHandler::as<std::vector<int16_t>>() const
    XDIAG_TRY_CATCH(return as_std_vector<int16_t>(key_, table_));

template <>
std::vector<int32_t> FileTomlHandler::as<std::vector<int32_t>>() const
    XDIAG_TRY_CATCH(return as_std_vector<int32_t>(key_, table_));

template <>
std::vector<int64_t> FileTomlHandler::as<std::vector<int64_t>>() const
    XDIAG_TRY_CATCH(return as_std_vector<int64_t>(key_, table_));

template <>
std::vector<uint8_t> FileTomlHandler::as<std::vector<uint8_t>>() const
    XDIAG_TRY_CATCH(return as_std_vector<uint8_t>(key_, table_));

template <>
std::vector<uint16_t> FileTomlHandler::as<std::vector<uint16_t>>() const
    XDIAG_TRY_CATCH(return as_std_vector<uint16_t>(key_, table_));

template <>
std::vector<uint32_t> FileTomlHandler::as<std::vector<uint32_t>>() const
    XDIAG_TRY_CATCH(return as_std_vector<uint32_t>(key_, table_));

template <>
std::vector<uint64_t> FileTomlHandler::as<std::vector<uint64_t>>() const
    XDIAG_TRY_CATCH(return as_std_vector<uint64_t>(key_, table_));

template <>
std::vector<double> FileTomlHandler::as<std::vector<double>>() const
    XDIAG_TRY_CATCH(return as_std_vector<double>(key_, table_));

template <>
std::vector<complex> FileTomlHandler::as<std::vector<complex>>() const
    XDIAG_TRY_CATCH(return as_std_vector<complex>(key_, table_));

template <>
std::vector<std::string> FileTomlHandler::as<std::vector<std::string>>() const
    XDIAG_TRY_CATCH(return as_std_vector<std::string>(key_, table_));

// Armadillo vectors
template <typename T>
static arma::Col<T> as_arma_vector(std::string key,
                                   toml::table const &table) try {
  auto node = table.at_path(key).node();
  if (node) {
    return arma_vector<T>(*node);
  } else {
    XDIAG_THROW(fmt::format("Key \"{}\" is not contained in TOML table", key));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
arma::vec FileTomlHandler::as<arma::vec>() const
    XDIAG_TRY_CATCH(return as_arma_vector<double>(key_, table_));

template <>
arma::cx_vec FileTomlHandler::as<arma::cx_vec>() const
    XDIAG_TRY_CATCH(return as_arma_vector<complex>(key_, table_));

template <>
arma::ivec FileTomlHandler::as<arma::ivec>() const
    XDIAG_TRY_CATCH(return as_arma_vector<arma::sword>(key_, table_));

template <>
arma::uvec FileTomlHandler::as<arma::uvec>() const
    XDIAG_TRY_CATCH(return as_arma_vector<arma::uword>(key_, table_));

// Armadillo matrices
template <typename T>
static arma::Mat<T> as_arma_matrix(std::string key,
                                   toml::table const &table) try {
  auto node = table.at_path(key).node();
  if (node) {
    return arma_matrix<T>(*node);
  } else {
    XDIAG_THROW(fmt::format("Key \"{}\" is not contained in TOML table", key));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
arma::mat FileTomlHandler::as<arma::mat>() const
    XDIAG_TRY_CATCH(return as_arma_matrix<double>(key_, table_));

template <>
arma::cx_mat FileTomlHandler::as<arma::cx_mat>() const
    XDIAG_TRY_CATCH(return as_arma_matrix<complex>(key_, table_));

template <>
arma::imat FileTomlHandler::as<arma::imat>() const
    XDIAG_TRY_CATCH(return as_arma_matrix<arma::sword>(key_, table_));

template <>
arma::umat FileTomlHandler::as<arma::umat>() const
    XDIAG_TRY_CATCH(return as_arma_matrix<arma::uword>(key_, table_));

template <> Permutation FileTomlHandler::as<Permutation>() const try {
  auto array = as_std_vector<int64_t>(key_, table_);
  return Permutation(array);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> PermutationGroup FileTomlHandler::as<PermutationGroup>() const try {
  auto node = table_.at_path(key_).node();
  if (node) {
    auto mat = arma_matrix<arma::sword>(*node);
    std::vector<Permutation> perms(mat.n_rows);
    for (std::size_t i = 0; i < mat.n_rows; ++i) {
      perms[i] =
          Permutation(std::vector<int64_t>(mat.begin_row(i), mat.end_row(i)));
    }
    return PermutationGroup(perms);
  } else {
    XDIAG_THROW(fmt::format("Key \"{}\" is not contained in TOML table", key_));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> Op FileTomlHandler::as<Op>() const try {
  auto node = table_.at_path(key_).node();
  if (node) {
    return op(*node);
  } else {
    XDIAG_THROW(fmt::format("Key \"{}\" is not contained in TOML table", key_));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> OpSum FileTomlHandler::as<OpSum>() const try {
  auto node = table_.at_path(key_).node();
  if (node) {
    return opsum(*node);
  } else {
    XDIAG_THROW(fmt::format("Key \"{}\" is not contained in TOML table", key_));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

//////////////////////////////////////////////////////////////////
// operator=

template <typename T> void FileTomlHandler::operator=(T const &value) try {
  table_.insert_or_assign(key_, value);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template void
    FileTomlHandler::operator= <std::string>(std::string const &value);
template void FileTomlHandler::operator= <int8_t>(int8_t const &value);
template void FileTomlHandler::operator= <int16_t>(int16_t const &value);
template void FileTomlHandler::operator= <int32_t>(int32_t const &value);
template void FileTomlHandler::operator= <int64_t>(int64_t const &value);
template void FileTomlHandler::operator= <uint8_t>(uint8_t const &value);
template void FileTomlHandler::operator= <uint16_t>(uint16_t const &value);
template void FileTomlHandler::operator= <uint32_t>(uint32_t const &value);
template void FileTomlHandler::operator= <double>(double const &value);

template <>
void FileTomlHandler::operator= <uint64_t>(uint64_t const &value) try {
  table_.insert_or_assign(key_, (int64_t)value);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> void FileTomlHandler::operator=(complex const &value) try {
  table_.insert_or_assign(key_, toml::array{value.real(), value.imag()});
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator=
    <std::vector<int8_t>>(std::vector<int8_t> const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<int16_t>>(std::vector<int16_t> const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<int32_t>>(std::vector<int32_t> const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<int64_t>>(std::vector<int64_t> const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint8_t>>(std::vector<uint8_t> const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint16_t>>(std::vector<uint16_t> const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint32_t>>(std::vector<uint32_t> const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint64_t>>(std::vector<uint64_t> const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<double>>(std::vector<double> const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<complex>>(std::vector<complex> const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator=
    <std::vector<std::string>>(std::vector<std::string> const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <arma::vec>(arma::vec const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <arma::cx_vec>(arma::cx_vec const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <arma::ivec>(arma::ivec const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <arma::uvec>(arma::uvec const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <arma::mat>(arma::mat const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator= <arma::cx_mat>(arma::cx_mat const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator= <arma::imat>(arma::imat const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <arma::umat>(arma::umat const &value) try {
  table_.insert_or_assign(key_, toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <Permutation>(Permutation const &value) try {
  table_.insert_or_assign(key_, toml_array(value.array()));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator=
    <PermutationGroup>(PermutationGroup const &group) try {
  arma::imat mat(group.size(), group.nsites());
  for (int64_t i = 0; i < group.size(); ++i) {
    auto int_vec = group[i].array();
    std::vector<arma::sword> vec(int_vec.begin(), int_vec.end());
    mat.row(i) = arma::irowvec(vec);
  }
  table_.insert_or_assign(key_, toml_array(mat));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> void FileTomlHandler::operator= <Op>(Op const &op) try {
  table_.insert_or_assign(key_, toml_array(op));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> void FileTomlHandler::operator= <OpSum>(OpSum const &ops) try {
  if (ops.constants().size() > 0) {
    table_.insert_or_assign(key_, toml_table(ops));
  } else {
    table_.insert_or_assign(key_, toml_array(ops));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::io
