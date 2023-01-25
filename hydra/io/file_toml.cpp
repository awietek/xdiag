#include "file_toml.h"

#include <filesystem>

#include <hydra/utils/logger.h>

namespace hydra {

FileToml::FileToml(std::string filename, char iomode)
    : FileToml(filename, std::string(1, iomode)) {}

FileToml::FileToml(std::string filename, std::string iomode)
    : filename_(filename), iomode_(iomode) {
  if ((iomode != "r") && (iomode != "w")) {
    Log.err("Error creating FileToml: iomode must be eiter \"r\" or \"w\"");
  }

  if (iomode == "r") {
    if (!std::filesystem::exists(filename)) {
      Log.err("Error reading FileToml: file {} for reading does not exist",
              filename);
    }

    try {
      table_ = toml::parse_file(filename);
    } catch (const toml::parse_error &err) {
      Log.err("Error reading FileToml: cannot parse {}: \n {}", filename,
              std::string(err.description()));
    }
  } else if (iomode == "w") {
    try {
      std::ofstream fl;
      fl.open(filename, std::ios::out);
      fl.close();
    } catch (std::exception const &e) {
      Log.err("Error creating FileToml: could not create {} for writing",
              filename);
    }
  }
}

FileToml::~FileToml() { write(); }

bool FileToml::defined(std::string key) const {
  return (bool)table_.at_path(key);
}

io::FileTomlHandler FileToml::operator[](std::string key) {
  return io::FileTomlHandler(key, table_);
}

void FileToml::write() const {
  if (iomode_ == "w") {
    std::ofstream fl;
    fl.open(filename_, std::ios::out);
    fl << table_;
    fl.close();
  }
}

void FileToml::close() const { write(); }

bool FileToml::operator==(FileToml const &other) const {
  return (filename_ == other.filename_) && (iomode_ == other.iomode_) &&
         (table_ == other.table_);
}

bool FileToml::operator!=(FileToml const &other) const {
  return !operator==(other);
}

} // namespace hydra
