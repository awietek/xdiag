#pragma once
#ifdef XDIAG_USE_HDF5

#include <string>

#include <hdf5.h>

namespace xdiag::hdf5 {

class FileH5Handler {
public:
  FileH5Handler() = delete;
  FileH5Handler(hid_t file_id, std::string field);
  FileH5Handler(FileH5Handler const &) = delete;
  FileH5Handler &operator=(FileH5Handler const &) = delete;

  template <class data_t> void operator=(data_t const &data);

private:
  hid_t file_id_;
  std::string field_;
};

} // namespace xdiag::hdf5

#endif
