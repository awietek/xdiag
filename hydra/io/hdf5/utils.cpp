#include "utils.h"
#ifdef HYDRA_USE_HDF5
#include <hydra/utils/logger.h>

namespace hydra::hdf5 {

std::string dataset_name(std::string full_name) {
  std::size_t loc;
  while ((loc = full_name.find("/")) != std::string::npos) {
    full_name = full_name.substr(loc + 1);
  }
  return full_name.empty() ? std::string("dataset") : full_name;
}

std::vector<hid_t> create_groups(hid_t file_id, std::string name,
                                 bool open_if_exists) {
  std::vector<hid_t> groups;
  std::size_t loc;
  while ((loc = name.find("/")) != std::string::npos) {
    // Create another group...

    // Ignore the first /, if there is a leading /.
    if (loc != 0) {
      hid_t pid = (groups.size() == 0) ? file_id : groups[groups.size() - 1];

      // Group already exists
      if (H5Lexists(pid, name.substr(0, loc).c_str(), H5P_DEFAULT) > 0) {
        hid_t gid = H5Gopen(pid, name.substr(0, loc).c_str(), H5P_DEFAULT);
        if (gid == H5I_INVALID_HID) {
          Log.err("Error in hydra hdf5: unable to open group \"{}\"",
                  name.substr(0, loc));
        }
        groups.push_back(gid);

        // Group doesn't exist -> it's created
      } else {
        hid_t gid = H5Gcreate(pid, name.substr(0, loc).c_str(), H5P_DEFAULT,
                              H5P_DEFAULT, H5P_DEFAULT);
        if (gid == H5I_INVALID_HID) {
          Log.err("Error in hydra hdf5: unable to  create group \"{}\"",
                  name.substr(0, loc));
        }
        groups.push_back(gid);
      }
    }
    name = name.substr(loc + 1);
  }
  return groups;
}

void close_groups(std::vector<hid_t> const &groups) {
  for (std::size_t i = 0; i < groups.size(); ++i) {
    H5Gclose(groups[i]);
  }
}

} // namespace hydra::hdf5
#endif
