#include "file_toml_handler.h"

#include <filesystem>

#include <hydra/utils/logger.h>

namespace hydra::io {

FileTomlHandler::FileTomlHandler(std::string key, FileToml &file)
    : key_(key), file_(file) {}

template <> std::string FileTomlHandler::as<std::string>() {
  return file_.table_[key_].as<std::string>;
}

template <> void FileTomlHandler::operator==(std::string value) {
  return file_.table_.insert_or_assign(key_, value);
}

} // namespace hydra::io
