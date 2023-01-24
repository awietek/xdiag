#pragma once

#include <hydra/io/file_toml.h>
#include <string>

namespace hydra::io {

class FileToml;
class FileTomlHandler {
public:
  FileTomlHandler() = default;
  FileTomlHandler(std::string key, FileToml & file);

  template <class data_t> data_t as();
  template <class data_t> void operator=(data_t const &data);

private:
  std::string key_;
  FileToml& file_;
};

} // namespace hydra::io
