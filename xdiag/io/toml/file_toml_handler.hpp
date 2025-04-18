// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

#include <xdiag/extern/toml++/toml.hpp>

namespace xdiag::io {

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

} // namespace xdiag::io
