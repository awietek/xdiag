// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_USE_HDF5

#include <string>
#include <vector>

#include <hdf5.h>

namespace xdiag::hdf5 {

std::string dataset_name(std::string full_name);
std::vector<hid_t> create_groups(hid_t file_id, std::string field);
void close_groups(std::vector<hid_t> const &groups);

} // namespace xdiag::hdf5
#endif
