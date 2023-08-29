#include <map>
#include <string>

#include <hydra/operators/bondlist.h>

namespace hydra::spinhalf {

const std::map<std::string, int64_t> special_bonds_nup = {
    {"HB", 0}, {"HEISENBERG", 0}, {"EXCHANGE", 0}, {"ISING", 0},
    {"SZ", 0}, {"S+", 1},         {"S-", -1},      {"SCALARCHIRALITY", 0}};

int64_t nup(BondList bonds, double precision = 1e-12);

} // namespace hydra::spinhalf
