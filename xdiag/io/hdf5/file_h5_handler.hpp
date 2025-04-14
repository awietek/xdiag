// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_USE_HDF5

#include <string>

#include <hdf5.h>
#include <xdiag/common.hpp>
#include <xdiag/io/hdf5/file_h5_subview.hpp>

namespace xdiag::hdf5 {

class FileH5Handler {
public:
  FileH5Handler() = delete;
  FileH5Handler(hid_t file_id, std::string field);
  FileH5Handler(FileH5Handler const &) = delete;
  FileH5Handler &operator=(FileH5Handler const &) = delete;

  template <class data_t> XDIAG_API void operator=(data_t const &data);

  hdf5::FileH5Submat col(int col_number);
  hdf5::FileH5Subcube slice(int slice_number);

private:
  hid_t file_id_;
  std::string field_;
};

} // namespace xdiag::hdf5

#endif
