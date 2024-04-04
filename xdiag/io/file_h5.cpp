#ifdef XDIAG_USE_HDF5

#include "file_h5.hpp"

#include <xdiag/utils/logger.hpp>

namespace xdiag {

FileH5::FileH5(std::string filename, char iomode)
    : FileH5(filename, std::string(1, iomode)) {}

FileH5::FileH5(std::string filename, std::string iomode)
    : filename_(filename), iomode_(iomode) {
  // Open file in read-only mode
  if (iomode == "r") {
    Log.err("Error in xdiag hdf5: read mode (r) currently not implemented");
  }

  // Open file in forced write mode
  else if (iomode == "w!") {
    file_id_ =
        H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id_ == H5I_INVALID_HID) {
      Log.err("Error in xdiag hdf5: can't open file in (forced) write mode: {}",
              filename);
    }
  }

  // Open file in secure write mode
  else if (iomode == "w") {
    file_id_ =
        H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id_ == H5I_INVALID_HID) {
      Log.err("Error in xdiag hdf5: can't open file in write mode: {}",
              filename);
    }
  }

  // Open file in append mode
  else if (iomode == "a") {
    Log.err("Error in xdiag hdf5: append mode (a) currently not implemented");

  } else {
    Log.err("Error in xdiag hdf5: invalid iomode, must be one of \"r\", \"w\", "
            "\"w!\", \"a\"");
  }
}

FileH5::~FileH5() { close(); }

hdf5::FileH5Handler FileH5::operator[](std::string key) {
  return hdf5::FileH5Handler(file_id_, key);
}

void FileH5::close() { H5Fclose(file_id_); }

bool FileH5::operator==(FileH5 const &other) const {
  return (filename_ == other.filename_) && (iomode_ == other.iomode_) &&
         (file_id_ == other.file_id_);
}

bool FileH5::operator!=(FileH5 const &other) const {
  return !operator==(other);
}

} // namespace xdiag

#endif
