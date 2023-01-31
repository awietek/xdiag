#pragma once

#include <string>
#include <vector>

#include <hdf5.h>

namespace hydra::hdf5 {

std::string dataset_name(std::string full_name);
std::vector<hid_t> create_groups(hid_t file_id, std::string field,
                                 bool open_if_exists = true);
void close_groups(std::vector<hid_t> const &groups);

} // namespace hydra::hdf5
