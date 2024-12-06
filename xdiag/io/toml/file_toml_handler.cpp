#include "file_toml_handler.hpp"

#include <array>
#include <complex>
#include <cstdint>
#include <vector>

#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/extern/toml++/toml.hpp>
#include <xdiag/io/toml/toml_conversion.hpp>
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
  auto node = table.at_path(key);
  if (node) {
    try {
      return get_toml_value<T>(node);
    } catch (Error const &e) {
      XDIAG_THROW(fmt::format("Key \"{}\" exists in TOML table but cannot be "
                              "converted to type \"{}\"",
                              key, utils::type_string<T>()));
    }
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
  auto node = table.at_path(key);
  if (node) {
    try {
      return toml_array_to_std_vector<T>(get_toml_array(node));
    } catch (Error const &e) {
      XDIAG_THROW(fmt::format("Key \"{}\" exists in TOML table but cannot be "
                              "converted to type \"{}\"",
                              key, utils::type_string<std::vector<T>>()));
    }
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
  auto node = table.at_path(key);
  if (node) {
    try {
      return toml_array_to_arma_vector<T>(get_toml_array(node));
    } catch (Error const &e) {
      XDIAG_THROW(fmt::format("Key \"{}\" exists in TOML table but cannot be "
                              "converted to type \"{}\"",
                              key, utils::type_string<arma::Col<T>>()));
    }
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
  auto node = table.at_path(key);
  if (node) {
    try {
      return toml_array_to_arma_matrix<T>(get_toml_array(node));
    } catch (Error const &e) {
      XDIAG_THROW(fmt::format("Key \"{}\" exists in TOML table but cannot be "
                              "converted to type \"{}\"",
                              key, utils::type_string<arma::Mat<T>>()));
    }
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
  auto array =
      toml_array_to_std_vector<int64_t>(get_toml_array(table_.at_path(key_)));
  return Permutation(array);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> PermutationGroup FileTomlHandler::as<PermutationGroup>() const try {
  auto mat = toml_array_to_arma_matrix<arma::sword>(
      get_toml_array(table_.at_path(key_)));
  std::vector<Permutation> perms(mat.n_rows);
  for (std::size_t i = 0; i < mat.n_rows; ++i) {
    perms[i] =
        Permutation(std::vector<int64_t>(mat.begin_row(i), mat.end_row(i)));
  }
  return PermutationGroup(perms);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> Representation FileTomlHandler::as<Representation>() const try {
  auto character_entry =
      table_.at_path(key_ + std::string(".characters")).as_array();
  if (character_entry) {
    auto characters = toml_array_to_std_vector<complex>(*character_entry);
    auto allowed_symmetries_entry =
        table_.at_path(key_ + std::string(".allowed_symmetries")).as_array();
    if (allowed_symmetries_entry) {
      auto allowed_symmetries =
          toml_array_to_std_vector<int64_t>(*allowed_symmetries_entry);
      return Representation(characters, allowed_symmetries);
    } else {
      return Representation(characters);
    }
  } else {
    XDIAG_THROW("Error reading Representation from toml file: no field "
                "\"characters\"!");
    return Representation();
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> Op FileTomlHandler::as<Op>() const try {
  return Op(toml_array_to_op(get_toml_array(table_.at_path(key_))));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> OpSum FileTomlHandler::as<OpSum>() const try {
  auto node = table_.at_path(key_);
  auto array_opt = node.as_array();
  auto table_opt = node.as_table();
  if (array_opt) {
    return toml_array_to_op_list(*array_opt);
  } else if (table_opt) {
    return toml_table_to_op_list(*table_opt);
  } else {
    Log.err("Error parsing toml file to OpSum: entry needs to be either an "
            "array or a table");
    return OpSum();
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
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<int16_t>>(std::vector<int16_t> const &value) try {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<int32_t>>(std::vector<int32_t> const &value) try {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<int64_t>>(std::vector<int64_t> const &value) try {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint8_t>>(std::vector<uint8_t> const &value) try {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint16_t>>(std::vector<uint16_t> const &value) try {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint32_t>>(std::vector<uint32_t> const &value) try {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint64_t>>(std::vector<uint64_t> const &value) try {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<double>>(std::vector<double> const &value) try {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator=
    <std::vector<complex>>(std::vector<complex> const &value) try {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator=
    <std::vector<std::string>>(std::vector<std::string> const &value) try {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <arma::vec>(arma::vec const &value) try {
  table_.insert_or_assign(key_, arma_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <arma::cx_vec>(arma::cx_vec const &value) try {
  table_.insert_or_assign(key_, arma_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <arma::ivec>(arma::ivec const &value) try {
  table_.insert_or_assign(key_, arma_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <arma::uvec>(arma::uvec const &value) try {
  table_.insert_or_assign(key_, arma_vector_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <arma::mat>(arma::mat const &value) try {
  table_.insert_or_assign(key_, arma_matrix_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator= <arma::cx_mat>(arma::cx_mat const &value) try {
  table_.insert_or_assign(key_, arma_matrix_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void FileTomlHandler::operator= <arma::imat>(arma::imat const &value) try {
  table_.insert_or_assign(key_, arma_matrix_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <arma::umat>(arma::umat const &value) try {
  table_.insert_or_assign(key_, arma_matrix_to_toml_array(value));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator= <Permutation>(Permutation const &value) try {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value.array()));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void FileTomlHandler::operator=
    <PermutationGroup>(PermutationGroup const &group) try {
  arma::imat mat(group.n_symmetries(), group.n_sites());
  for (int64_t i = 0; i < group.n_symmetries(); ++i) {
    auto int_vec = group[i].array();
    std::vector<arma::sword> vec(int_vec.begin(), int_vec.end());
    mat.row(i) = arma::irowvec(vec);
  }
  table_.insert_or_assign(key_, arma_matrix_to_toml_array(mat));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void
    FileTomlHandler::operator= <Representation>(Representation const &rep) try {
  toml::table rep_table;
  rep_table.insert_or_assign("characters",
                             std_vector_to_toml_array(rep.characters()));
  rep_table.insert_or_assign(
      "allowed_symmetries", std_vector_to_toml_array(rep.allowed_symmetries()));
  table_.insert_or_assign(key_, rep_table);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> void FileTomlHandler::operator= <Op>(Op const &op) try {
  table_.insert_or_assign(key_, op_to_toml_array(op));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> void FileTomlHandler::operator= <OpSum>(OpSum const &ops) try {
  if (ops.constants().size() > 0) {
    table_.insert_or_assign(key_, op_list_to_toml_table(ops));
  } else {
    table_.insert_or_assign(key_, op_list_to_toml_array(ops));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::io
