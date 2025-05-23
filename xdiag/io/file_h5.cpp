// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#ifdef XDIAG_USE_HDF5

#include "file_h5.hpp"

#include <xdiag/utils/logger.hpp>

namespace xdiag {

FileH5::FileH5(std::string filename, char iomode)
    : FileH5(filename, std::string(1, iomode)) {}

FileH5::FileH5(std::string filename, std::string iomode) try
    : filename_(filename), iomode_(iomode), closed_(false) {
  // Open file in read-only mode
  if (iomode == "r") {
    Log(2, "opening h5file in r mode.");
    file_id_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    XDIAG_THROW("Error in xdiag hdf5: read mode (r) currently not implemented");
    if (file_id_ == H5I_INVALID_HID) {
      XDIAG_THROW(fmt::format(
          "Cannot open file in read mode \"r\": {}\n Maybe it does not exist?",
          filename));
    }
  }

  // Open file in forced write mode
  else if (iomode == "w!") {
    Log(2, "creating h5file in w! mode.");
    file_id_ =
        H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id_ == H5I_INVALID_HID) {
      XDIAG_THROW(fmt::format(
          "Cannot open file in forced write mode \"w!\": {}", filename));
    }
  }

  // Open file in secure write mode
  else if (iomode == "w") {
    Log(2, "creating h5file in w mode.");
    file_id_ =
        H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id_ == H5I_INVALID_HID) {
      XDIAG_THROW(fmt::format("Cannot open file in secure write mode \"w\": "
                              "{}\n Maybe it already exists?",
                              filename));
    }
  }

  // Open file in append mode
  else if (iomode == "a") {
    Log(2, "opening h5file in append mode.");
    file_id_ = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id_ == H5I_INVALID_HID) {
      XDIAG_THROW(
          fmt::format("Cannot open file in append mode \"w!\": {}", filename));
    }
  } else {
    XDIAG_THROW(
        "Error in xdiag hdf5: invalid iomode, must be one of \"r\", \"w\", "
        "\"w!\", \"a\"");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

FileH5::~FileH5() {
  if (!closed_) {
    close();
  }
  closed_ = true;
}

void FileH5::close() {
  H5Fclose(file_id_);
  closed_ = true;
}

hdf5::FileH5Handler FileH5::operator[](std::string key) {
  return hdf5::FileH5Handler(file_id_, key);
}

bool FileH5::operator==(FileH5 const &other) const {
  return (filename_ == other.filename_) && (iomode_ == other.iomode_) &&
         (file_id_ == other.file_id_);
}

bool FileH5::operator!=(FileH5 const &other) const {
  return !operator==(other);
}

} // namespace xdiag

#endif
