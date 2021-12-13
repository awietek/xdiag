#pragma once

#include <vector>

namespace hydra::symmetries {

// Creates the table of Fermi signs for a given particle number
template <typename bit_t, class GroupAction>
std::vector<bool> fermi_bool_table(int npar, GroupAction const &group_action);

} // namespace hydra::symmetries
