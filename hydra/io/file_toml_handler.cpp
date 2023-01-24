#include "file_toml_handler.h"

#include <array>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <vector>

#include <hydra/utils/logger.h>

namespace hydra::io {

FileTomlHandler::FileTomlHandler(std::string key, toml::table &table)
    : key_(key), table_(table) {}

template <typename T>
static T retrieve_plain(std::string const &key, toml::table const &table) {
  auto val = table.at_path(key).value<T>();
  if (val) {
    return T(*val);
  } else {
    Log.err("Error parsing toml: key \"{}\" not found!", key);
    return T();
  }
}

template <typename T, int size>
static std::array<T, size> retrieve_std_array(std::string const &key,
                                              toml::table const &table) {
  auto arr = table.at_path(key).as_array();
  if (arr) {
    int i = 0;
    std::array<T, size> array;
    for (auto &&x : *arr) {
      if (i < size) {
        array[i] = *(x.value<T>());
      }
      ++i;
    }
    if (i != size) {
      Log.err(
          "Error reading TOML node: array size ({}) exceeds expectation ({})!");
    }
    return array;
  } else {
    Log.err("Error parsing toml: key \"{}\" not found!", key);
    return std::array<T, size>();
  }
}

template <typename T>
static std::vector<T> retrieve_std_vector(std::string const &key,
                                          toml::table const &table) {
  auto arr = table.at_path(key).as_array();
  if (arr) {
    std::vector<T> vector(arr->size());
    int i = 0;
    for (auto &&x : *arr) {
      vector[i] = *(x.value<T>());
      ++i;
    }
    return vector;
  } else {
    Log.err("Error parsing toml: key \"{}\" not found!", key);
    return std::vector<T>();
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

template <> void FileTomlHandler::operator=(std::string const &value) {
  table_.insert_or_assign(key_, value);
}

} // namespace hydra::io
