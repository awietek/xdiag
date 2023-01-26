#pragma once

#include <extern/toml++/toml.h>

#include <string>

namespace hydra::io {

class FileTomlHandler {
public:
  FileTomlHandler(std::string key, toml::table &file);
  FileTomlHandler(FileTomlHandler const &) = delete;
  FileTomlHandler &operator=(FileTomlHandler const &) = delete;

  template <class data_t> data_t as() const;
  template <class data_t> void operator=(data_t const &data);

private:
  std::string key_;
  toml::table &table_;
};

} // namespace hydra::io
