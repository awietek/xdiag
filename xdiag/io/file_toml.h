#pragma once

#include <xdiag/extern/toml++/toml.h>
#include <xdiag/io/toml/file_toml_handler.h>
#include <string>

namespace xdiag {

class FileToml {
public:
  FileToml() = default;
  ~FileToml();
  FileToml(std::string filename, std::string iomode = "r");
  FileToml(std::string filename, char iomode);

  bool defined(std::string key) const;
  void write() const;
  void close() const;

  io::FileTomlHandler operator[](std::string key);
  bool operator==(FileToml const &other) const;
  bool operator!=(FileToml const &other) const;

  toml::table table() const;
  
private:
  std::string filename_;
  std::string iomode_;
  toml::table table_;
};

} // namespace xdiag
