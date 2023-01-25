#pragma once

#include <extern/toml++/toml.h>
#include <hydra/io/file_toml_handler.h>
#include <string>

namespace hydra {

class FileToml {
public:
  FileToml() = default;
  ~FileToml();
  FileToml(std::string filename, char iomode = 'r');
  FileToml(std::string filename, std::string iomode = "r");

  bool defined(std::string key) const;
  io::FileTomlHandler operator[](std::string key);
  void write() const;
  void close() const;

  bool operator==(FileToml const &other) const;
  bool operator!=(FileToml const &other) const;

private:
  std::string filename_;
  std::string iomode_;
  toml::table table_;
};

} // namespace hydra
