#include <map>
#include <string>

namespace hydra::terms::spinhalf {

const std::map<std::string, int> special_bonds_nup = {
    {"HB", 0}, {"HEISENBERG", 0}, {"EXCHANGE", 0}, {"ISING", 0},
    {"SZ", 0}, {"S+", 1},         {"S-", -1},      {"SCALARCHIRALITY"}};

  
} // namespace hydra::terms::spinhalf
