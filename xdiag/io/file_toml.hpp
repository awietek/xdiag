#pragma once

#include <istream>
#include <string>

#include <xdiag/common.hpp>

#include <xdiag/extern/toml++/toml.hpp>
#include <xdiag/io/toml/file_toml_handler.hpp>

namespace xdiag {

class FileToml {
public:
  XDIAG_API FileToml() = default;
  XDIAG_API FileToml(const char *filename);
  XDIAG_API FileToml(std::string filename);
  XDIAG_API FileToml(std::istream &is);

  XDIAG_API bool defined(std::string key) const;
  XDIAG_API io::FileTomlHandler operator[](std::string key);
  XDIAG_API void write(std::string filename, std::string mode = "w") const;

  XDIAG_API bool operator==(FileToml const &other) const;
  XDIAG_API bool operator!=(FileToml const &other) const;

  std::vector<std::string> keys() const;
  toml::table table() const;

private:
  toml::table table_;
};

} // namespace xdiag
