// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_USE_HDF5

#include <hdf5.h>
#include <map>
#include <string>
#include <vector>

#include <xdiag/io/hdf5/file_h5_handler.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

class XDIAG_API FileH5 {
public:
  FileH5() = default;
  FileH5(std::string filename, char iomode);
  FileH5(std::string filename, std::string iomode = "w");
  ~FileH5();

  void close();

  hdf5::FileH5Handler operator[](std::string key);
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
