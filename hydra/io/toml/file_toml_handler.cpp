#include "file_toml_handler.h"

#include <extern/armadillo/armadillo>
#include <extern/toml++/toml.h>

#include <array>
#include <complex>
#include <cstdint>
#include <vector>

#include <hydra/common.h>
#include <hydra/io/toml/toml_conversion.h>
#include <hydra/utils/logger.h>

#include <hydra/symmetries/permutation.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

#include <hydra/operators/bond.h>
#include <hydra/operators/bondlist.h>

namespace hydra::io {

FileTomlHandler::FileTomlHandler(std::string key, toml::table &table)
    : key_(key), table_(table) {}

// Plain values
template <> bool FileTomlHandler::as<bool>() const {
  return get_toml_value<bool>(table_.at_path(key_));
}

template <> int8_t FileTomlHandler::as<int8_t>() const {
  return get_toml_value<int8_t>(table_.at_path(key_));
}

template <> int16_t FileTomlHandler::as<int16_t>() const {
  return get_toml_value<int16_t>(table_.at_path(key_));
}

template <> int32_t FileTomlHandler::as<int32_t>() const {
  return get_toml_value<int32_t>(table_.at_path(key_));
}

template <> int64_t FileTomlHandler::as<int64_t>() const {
  return get_toml_value<int64_t>(table_.at_path(key_));
}

template <> uint8_t FileTomlHandler::as<uint8_t>() const {
  return get_toml_value<uint8_t>(table_.at_path(key_));
}

template <> uint16_t FileTomlHandler::as<uint16_t>() const {
  return get_toml_value<uint16_t>(table_.at_path(key_));
}

template <> uint32_t FileTomlHandler::as<uint32_t>() const {
  return get_toml_value<uint32_t>(table_.at_path(key_));
}

template <> uint64_t FileTomlHandler::as<uint64_t>() const {
  return get_toml_value<uint64_t>(table_.at_path(key_));
}

template <> double FileTomlHandler::as<double>() const {
  return get_toml_value<double>(table_.at_path(key_));
}

template <> complex FileTomlHandler::as<complex>() const {
  return get_toml_value<complex>(table_.at_path(key_));
}

template <> std::string FileTomlHandler::as<std::string>() const {
  return get_toml_value<std::string>(table_.at_path(key_));
}

// std::vectors
template <>
std::vector<int8_t> FileTomlHandler::as<std::vector<int8_t>>() const {
  return toml_array_to_std_vector<int8_t>(get_toml_array(table_.at_path(key_)));
}
template <>
std::vector<int16_t> FileTomlHandler::as<std::vector<int16_t>>() const {
  return toml_array_to_std_vector<int16_t>(
      get_toml_array(table_.at_path(key_)));
}
template <>
std::vector<int32_t> FileTomlHandler::as<std::vector<int32_t>>() const {
  return toml_array_to_std_vector<int32_t>(
      get_toml_array(table_.at_path(key_)));
}
template <>
std::vector<int64_t> FileTomlHandler::as<std::vector<int64_t>>() const {
  return toml_array_to_std_vector<int64_t>(
      get_toml_array(table_.at_path(key_)));
}
template <>
std::vector<uint8_t> FileTomlHandler::as<std::vector<uint8_t>>() const {
  return toml_array_to_std_vector<uint8_t>(
      get_toml_array(table_.at_path(key_)));
}
template <>
std::vector<uint16_t> FileTomlHandler::as<std::vector<uint16_t>>() const {
  return toml_array_to_std_vector<uint16_t>(
      get_toml_array(table_.at_path(key_)));
}
template <>
std::vector<uint32_t> FileTomlHandler::as<std::vector<uint32_t>>() const {
  return toml_array_to_std_vector<uint32_t>(
      get_toml_array(table_.at_path(key_)));
}
template <>
std::vector<uint64_t> FileTomlHandler::as<std::vector<uint64_t>>() const {
  return toml_array_to_std_vector<uint64_t>(
      get_toml_array(table_.at_path(key_)));
}

template <>
std::vector<double> FileTomlHandler::as<std::vector<double>>() const {
  return toml_array_to_std_vector<double>(get_toml_array(table_.at_path(key_)));
}
template <>
std::vector<complex> FileTomlHandler::as<std::vector<complex>>() const {
  return toml_array_to_std_vector<complex>(
      get_toml_array(table_.at_path(key_)));
}
template <>
std::vector<std::string> FileTomlHandler::as<std::vector<std::string>>() const {
  return toml_array_to_std_vector<std::string>(
      get_toml_array(table_.at_path(key_)));
}

// Armadillo vectors
template <> arma::vec FileTomlHandler::as<arma::vec>() const {
  return toml_array_to_arma_vector<double>(
      get_toml_array(table_.at_path(key_)));
}

template <> arma::cx_vec FileTomlHandler::as<arma::cx_vec>() const {
  return toml_array_to_arma_vector<complex>(
      get_toml_array(table_.at_path(key_)));
}

template <> arma::ivec FileTomlHandler::as<arma::ivec>() const {
  return toml_array_to_arma_vector<arma::sword>(
      get_toml_array(table_.at_path(key_)));
}

template <> arma::uvec FileTomlHandler::as<arma::uvec>() const {
  return toml_array_to_arma_vector<arma::uword>(
      get_toml_array(table_.at_path(key_)));
}

// Armadillo matrices
template <> arma::mat FileTomlHandler::as<arma::mat>() const {
  return toml_array_to_arma_matrix<double>(
      get_toml_array(table_.at_path(key_)));
}

template <> arma::cx_mat FileTomlHandler::as<arma::cx_mat>() const {
  return toml_array_to_arma_matrix<complex>(
      get_toml_array(table_.at_path(key_)));
}

template <> arma::imat FileTomlHandler::as<arma::imat>() const {
  return toml_array_to_arma_matrix<arma::sword>(
      get_toml_array(table_.at_path(key_)));
}

template <> arma::umat FileTomlHandler::as<arma::umat>() const {
  return toml_array_to_arma_matrix<arma::uword>(
      get_toml_array(table_.at_path(key_)));
}

template <> Permutation FileTomlHandler::as<Permutation>() const {
  auto array =
      toml_array_to_std_vector<int64_t>(get_toml_array(table_.at_path(key_)));
  return Permutation(array);
}

template <> PermutationGroup FileTomlHandler::as<PermutationGroup>() const {
  auto mat = toml_array_to_arma_matrix<arma::sword>(
      get_toml_array(table_.at_path(key_)));
  std::vector<Permutation> perms(mat.n_rows);
  for (std::size_t i = 0; i < mat.n_rows; ++i) {
    perms[i] = Permutation(std::vector<int64_t>(mat.begin_row(i), mat.end_row(i)));
  }
  return PermutationGroup(perms);
}

template <> Representation FileTomlHandler::as<Representation>() const {
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
    Log.err("Error reading Representation from toml file: no field "
            "\"characters\"!");
    return Representation();
  }
}

template <> Bond FileTomlHandler::as<Bond>() const {
  return Bond(toml_array_to_bond(get_toml_array(table_.at_path(key_))));
}

template <> BondList FileTomlHandler::as<BondList>() const {
  auto node = table_.at_path(key_);
  auto array_opt = node.as_array();
  auto table_opt = node.as_table();
  if (array_opt) {
    return toml_array_to_bond_list(*array_opt);
  } else if (table_opt) {
    return toml_table_to_bond_list(*table_opt);
  } else {
    Log.err("Error parsing toml file to BondList: entry needs to be either an "
            "array or a table");
    return BondList();
  }
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
  table_.insert_or_assign(key_, toml::array{value.real(), value.imag()});
}

template <>
void FileTomlHandler::operator=
    <std::vector<int8_t>>(std::vector<int8_t> const &value) {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
}
template <>
void FileTomlHandler::operator=
    <std::vector<int16_t>>(std::vector<int16_t> const &value) {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
}
template <>
void FileTomlHandler::operator=
    <std::vector<int32_t>>(std::vector<int32_t> const &value) {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
}
template <>
void FileTomlHandler::operator=
    <std::vector<int64_t>>(std::vector<int64_t> const &value) {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint8_t>>(std::vector<uint8_t> const &value) {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint16_t>>(std::vector<uint16_t> const &value) {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint32_t>>(std::vector<uint32_t> const &value) {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint64_t>>(std::vector<uint64_t> const &value) {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
}
template <>
void FileTomlHandler::operator=
    <std::vector<double>>(std::vector<double> const &value) {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
}
template <>
void FileTomlHandler::operator=
    <std::vector<complex>>(std::vector<complex> const &value) {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
}

template <>
void FileTomlHandler::operator=
    <std::vector<std::string>>(std::vector<std::string> const &value) {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value));
}

template <> void FileTomlHandler::operator=<arma::vec>(arma::vec const &value) {
  table_.insert_or_assign(key_, arma_vector_to_toml_array(value));
}

template <>
void FileTomlHandler::operator=<arma::cx_vec>(arma::cx_vec const &value) {
  table_.insert_or_assign(key_, arma_vector_to_toml_array(value));
}

template <>
void FileTomlHandler::operator=<arma::ivec>(arma::ivec const &value) {
  table_.insert_or_assign(key_, arma_vector_to_toml_array(value));
}

template <>
void FileTomlHandler::operator=<arma::uvec>(arma::uvec const &value) {
  table_.insert_or_assign(key_, arma_vector_to_toml_array(value));
}

template <> void FileTomlHandler::operator=<arma::mat>(arma::mat const &value) {
  table_.insert_or_assign(key_, arma_matrix_to_toml_array(value));
}
template <>
void FileTomlHandler::operator=<arma::cx_mat>(arma::cx_mat const &value) {
  table_.insert_or_assign(key_, arma_matrix_to_toml_array(value));
}
template <>
void FileTomlHandler::operator=<arma::imat>(arma::imat const &value) {
  table_.insert_or_assign(key_, arma_matrix_to_toml_array(value));
}

template <>
void FileTomlHandler::operator=<arma::umat>(arma::umat const &value) {
  table_.insert_or_assign(key_, arma_matrix_to_toml_array(value));
}

template <>
void FileTomlHandler::operator=<Permutation>(Permutation const &value) {
  table_.insert_or_assign(key_, std_vector_to_toml_array(value.array()));
}

template <>
void FileTomlHandler::operator=
    <PermutationGroup>(PermutationGroup const &group) {
  arma::imat mat(group.n_symmetries(), group.n_sites());
  for (int64_t i = 0; i < group.n_symmetries(); ++i) {
    auto int_vec = group[i].array();
    std::vector<arma::sword> vec(int_vec.begin(), int_vec.end());
    mat.row(i) = arma::irowvec(vec);
  }
  table_.insert_or_assign(key_, arma_matrix_to_toml_array(mat));
}

template <>
void FileTomlHandler::operator=<Representation>(Representation const &rep) {
  toml::table rep_table;
  rep_table.insert_or_assign("characters",
                             std_vector_to_toml_array(rep.characters()));
  rep_table.insert_or_assign(
      "allowed_symmetries", std_vector_to_toml_array(rep.allowed_symmetries()));
  table_.insert_or_assign(key_, rep_table);
}

template <> void FileTomlHandler::operator=<Bond>(Bond const &bond) {
  table_.insert_or_assign(key_, bond_to_toml_array(bond));
}

template <> void FileTomlHandler::operator=<BondList>(BondList const &bonds) {
  if ((bonds.couplings().size() > 0) || (bonds.matrices().size() > 0)) {
    table_.insert_or_assign(key_, bond_list_to_toml_table(bonds));
  } else {
    table_.insert_or_assign(key_, bond_list_to_toml_array(bonds));
  }
}

} // namespace hydra::io
