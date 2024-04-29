#pragma once
#ifdef XDIAG_USE_HDF5

#include <string>

#include <hdf5.h>

namespace xdiag::hdf5 {

class FileH5Submat {
public:
  FileH5Submat() = delete;
  FileH5Submat(hid_t file_id, std::string field, int row_number, int col_number);
  FileH5Submat(FileH5Submat const &) = delete;
  FileH5Submat &operator=(FileH5Submat const &) = delete;
  template <class data_t> void operator=(data_t const &data);

private:
  hid_t file_id_;
  std::string field_;
  int col_number_;
  int row_number_;
};


class FileH5Subcube {
public:
  FileH5Subcube() = delete;
  FileH5Subcube(hid_t file_id, std::string field, int row_number, int col_number, int slice_number);
  FileH5Subcube(FileH5Subcube const &) = delete;
  FileH5Subcube &operator=(FileH5Subcube const &) = delete;
  template <class data_t> void operator=(data_t const &data);

private:
  hid_t file_id_;
  std::string field_;
  int col_number_;
  int row_number_;
  int slice_number_;
};


} // namespace xdiag::hdf5

#endif
