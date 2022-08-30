#pragma once

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/operators/operator_utils.h>

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace hydra::utils {

const std::vector<std::string> spinhalf_bonds = {
    "HB", "HEISENBERG", "EXCHANGE", "ISING", "SZ", "SX", "SY", "S+", "S-"};
const std::map<std::string, int> spinhalf_bond_nup = {
    {"HB", 0},  {"HEISENBERG", 0},    {"EXCHANGE", 0},      {"ISING", 0},
    {"SZ", 0},  {"SX", undefined_qn}, {"SY", undefined_qn}, {"S+", 1},
    {"S-", -1},
};

const std::vector<std::string> electron_bonds = {
    "HB", "HEISENBERG", "EXCHANGE", "ISING", "SZ",     "SX", "SY",
    "S+", "S-",         "Cdagup",   "Cup",   "Cdagdn", "Cdn"};
const std::map<std::string, std::pair<int, int>> electron_bond_nup_ndn = {
    {"HB", {0, 0}},        {"HEISENBERG", {0, 0}}, {"EXCHANGE", {0, 0}},
    {"ISING", {0, 0}},     {"SZ", {0, 0}},         {"SX", undefined_qns},
    {"SY", undefined_qns}, {"S+", {1, -1}},        {"S-", {-1, 1}},
    {"Cdagup", {1, 0}},    {"Cup", {-1, 0}},       {"Cdagdn", {0, 1}},
    {"Cdn", {0, -1}}};

int spinhalf_nup(BondList const &bonds, Couplings const &cpls);
std::pair<int, int> tj_nup_ndn(BondList const &bonds, Couplings const &cpls);
std::pair<int, int> electron_nup_ndn(BondList const &bonds,
                                     Couplings const &cpls);

template <class Block>
inline int spinhalf_nup(BondList const &bonds, Couplings const &cpls,
                        Block &&spinhalf) {
  int op_nup = spinhalf_nup(bonds, cpls);
  if ((op_nup == undefined_qn) || (spinhalf.n_up() == undefined_qn)) {
    return undefined_qn;
  } else {
    return op_nup + spinhalf.n_up();
  }
}

template <class Block>
inline std::pair<int, int> tj_nup_ndn(BondList const &bonds,
                                      Couplings const &cpls, Block &&tj) {
  auto op_ns = tj_nup_ndn(bonds, cpls);
  if ((op_ns == undefined_qns) || (tj.n_up() == undefined_qn) ||
      (tj.n_dn() == undefined_qn)) {
    return undefined_qns;
  } else {
    return {op_ns.first + tj.n_up(), op_ns.second + tj.n_dn()};
  }
}
template <class Block>
inline std::pair<int, int> electron_nup_ndn(BondList const &bonds,
                                            Couplings const &cpls,
                                            Block &&electron) {
  auto op_ns = electron_nup_ndn(bonds, cpls);
  if ((op_ns == undefined_qns) || (electron.n_up() == undefined_qn) ||
      (electron.n_dn() == undefined_qn)) {
    return undefined_qns;
  } else {
    return {op_ns.first + electron.n_up(), op_ns.second + electron.n_dn()};
  }
}

} // namespace hydra::utils
