#include "file_toml.hpp"

#include <fstream>
#include <xdiag/extern/fmt/format.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

static toml::table istream_to_toml_table(std::istream &is) try {
  if (!is.good()) {
    XDIAG_THROW("Unable to read FileToml. Invalid input stream.");
  }
  try {
    return toml::parse(is);
  } catch (toml::parse_error const &err) {
    XDIAG_THROW(fmt::format("Unable to read toml input stream:\n{}",
                            err.description()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

FileToml::FileToml(const char *filename) try {
  std::ifstream f(filename);
  if (!f.good()) {
    XDIAG_THROW(fmt::format(
        "Unable to read FileToml. File \"{}\" does not exist", filename));
  }
  table_ = istream_to_toml_table(f);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

FileToml::FileToml(std::string filename) try : FileToml(filename.c_str()) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

FileToml::FileToml(std::istream &is) try {
  if (!is.good()) {
    XDIAG_THROW("Unable to read FileToml. Invalid input stream.");
  }
  table_ = istream_to_toml_table(is);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool FileToml::defined(std::string key) const try {
  return (bool)table_.at_path(key);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

io::FileTomlHandler FileToml::operator[](std::string key) try {
  return io::FileTomlHandler(key, table_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void FileToml::write(std::string filename, std::string mode) const try {
  std::ofstream fl;
  if (mode == "w") {
    fl.open(filename, std::ios::out);
  } else if (mode == "a") {
    fl.open(filename, std::ios::app);
  } else {
    XDIAG_THROW("Invalid mode for writing. Must be one of \"w\" or \"a\"");
  }
  if (fl.is_open()) {
    fl << table_;
    fl.close();
  } else {
    XDIAG_THROW(fmt::format("Unable to open file \"{}\"", filename));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool FileToml::operator==(FileToml const &other) const {
  return table_ == other.table_;
}

bool FileToml::operator!=(FileToml const &other) const {
  return !operator==(other);
}

std::vector<std::string> FileToml::keys() const {
  std::vector<std::string> ks;
  for (auto const &kv : table_) {
    ks.push_back(std::string(kv.first.str()));
  }
  return ks;
}
toml::table FileToml::table() const { return table_; }

} // namespace xdiag
