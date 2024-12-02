#pragma once

#include <xdiag/extern/toml++/toml.hpp>
#include <xdiag/io/toml/file_toml_handler.hpp>

#include <istream>
#include <string>

namespace xdiag {

class FileToml {
public:
  FileToml() = default;
  FileToml(const char* filename);
  FileToml(std::string filename);
  FileToml(std::istream &is);

  bool defined(std::string key) const;
  io::FileTomlHandler operator[](std::string key);

  void write(std::string filename, std::string mode = "w") const;

  bool operator==(FileToml const &other) const;
  bool operator!=(FileToml const &other) const;

  std::vector<std::string> keys() const;
  toml::table table() const;

private:
  toml::table table_;
};

} // namespace xdiag
