#include "file_toml_handler.h"

#include <extern/armadillo/armadillo>
#include <extern/toml++/toml.h>

#include <array>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <vector>

#include <hydra/common.h>
#include <hydra/io/toml_conversion.h>
#include <hydra/utils/logger.h>

#include <hydra/symmetries/permutation.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

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
      toml_array_to_std_vector<int>(get_toml_array(table_.at_path(key_)));
  return Permutation(array);
}

//////////////////////////////////////////////////////////////////
// operator=

template <typename T> void FileTomlHandler::operator=(T const &value) {
  insert(key_, value, table_);
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
  insert(key_, (int64_t)value, table_);
}

template <> void FileTomlHandler::operator=(complex const &value) {
  insert(key_, toml::array{value.real(), value.imag()}, table_);
}

template <>
void FileTomlHandler::operator=
    <std::vector<int8_t>>(std::vector<int8_t> const &value) {
  insert_std_vector(key_, value, table_);
}
template <>
void FileTomlHandler::operator=
    <std::vector<int16_t>>(std::vector<int16_t> const &value) {
  insert_std_vector(key_, value, table_);
}
template <>
void FileTomlHandler::operator=
    <std::vector<int32_t>>(std::vector<int32_t> const &value) {
  insert_std_vector(key_, value, table_);
}
template <>
void FileTomlHandler::operator=
    <std::vector<int64_t>>(std::vector<int64_t> const &value) {
  insert_std_vector(key_, value, table_);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint8_t>>(std::vector<uint8_t> const &value) {
  insert_std_vector(key_, value, table_);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint16_t>>(std::vector<uint16_t> const &value) {
  insert_std_vector(key_, value, table_);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint32_t>>(std::vector<uint32_t> const &value) {
  insert_std_vector(key_, value, table_);
}
template <>
void FileTomlHandler::operator=
    <std::vector<uint64_t>>(std::vector<uint64_t> const &value) {
  insert_std_vector(key_, value, table_);
}
template <>
void FileTomlHandler::operator=
    <std::vector<double>>(std::vector<double> const &value) {
  insert_std_vector(key_, value, table_);
}
template <>
void FileTomlHandler::operator=
    <std::vector<complex>>(std::vector<complex> const &value) {
  insert_std_vector(key_, value, table_);
}

template <>
void FileTomlHandler::operator=
    <std::vector<std::string>>(std::vector<std::string> const &value) {
  insert_std_vector(key_, value, table_);
}

template <> void FileTomlHandler::operator=<arma::vec>(arma::vec const &value) {
  insert_arma_vector(key_, value, table_);
}

template <>
void FileTomlHandler::operator=<arma::cx_vec>(arma::cx_vec const &value) {
  insert_arma_vector(key_, value, table_);
}

template <>
void FileTomlHandler::operator=<arma::ivec>(arma::ivec const &value) {
  insert_arma_vector(key_, value, table_);
}

template <>
void FileTomlHandler::operator=<arma::uvec>(arma::uvec const &value) {
  insert_arma_vector(key_, value, table_);
}

template <> void FileTomlHandler::operator=<arma::mat>(arma::mat const &value) {
  insert_arma_matrix(key_, value, table_);
}
template <>
void FileTomlHandler::operator=<arma::cx_mat>(arma::cx_mat const &value) {
  insert_arma_matrix(key_, value, table_);
}
template <>
void FileTomlHandler::operator=<arma::imat>(arma::imat const &value) {
  insert_arma_matrix(key_, value, table_);
}

template <>
void FileTomlHandler::operator=<arma::umat>(arma::umat const &value) {
  insert_arma_matrix(key_, value, table_);
}

template <>
void FileTomlHandler::operator=<Permutation>(Permutation const &value) {
  insert_std_vector(key_, value.array(), table_);
}
  
} // namespace hydra::io
