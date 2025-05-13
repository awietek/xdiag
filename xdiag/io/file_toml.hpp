// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

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
  XDIAG_API explicit FileToml(const char *filename);
  XDIAG_API explicit FileToml(std::string filename);
  XDIAG_API explicit FileToml(std::istream &is);

  bool defined(std::string key) const;
  XDIAG_API io::FileTomlHandler operator[](std::string key);
  void write(std::string filename, std::string mode = "w") const;

  XDIAG_API bool operator==(FileToml const &other) const;
  XDIAG_API bool operator!=(FileToml const &other) const;

  std::vector<std::string> keys() const;
  toml::table table() const;

private:
  toml::table table_;
};

XDIAG_API bool defined(FileToml const &fl, std::string key);

} // namespace xdiag
