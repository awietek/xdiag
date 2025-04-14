// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_USE_HDF5

#include <map>
#include <string>
#include <vector>

#include <hdf5.h>

#include <xdiag/common.hpp>
#include <xdiag/io/hdf5/file_h5_handler.hpp>

namespace xdiag {

class FileH5 {
public:
  FileH5() = default;
  XDIAG_API FileH5(std::string filename, char iomode);
  XDIAG_API FileH5(std::string filename, std::string iomode = "w");
  XDIAG_API ~FileH5();

  void close();

  XDIAG_API hdf5::FileH5Handler operator[](std::string key);
  bool operator==(FileH5 const &other) const;
  bool operator!=(FileH5 const &other) const;

private:
  std::string filename_;
  std::string iomode_;
  hid_t file_id_;
  bool closed_;
};

} // namespace xdiag

#endif
