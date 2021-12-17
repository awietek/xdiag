#pragma once

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/operators/operator_utils.h>

#include <limits>
#include <map>
#include <string>
#include <utility>

namespace hydra::utils {

constexpr int undefined_qn = std::numeric_limits<int>::min();
constexpr std::pair<int, int> undefined_qns = {undefined_qn, undefined_qn};

const std::map<std::string, int> spinhalf_bond_nup = {
    {"HB", 0},  {"HEISENBERG", 0},    {"EXCHANGE", 0},      {"ISING", 0},
    {"SZ", 0},  {"SX", undefined_qn}, {"SY", undefined_qn}, {"S+", 1},
    {"S-", -1},
};

int spinhalf_operator_nup(BondList const &bonds, Couplings const &cpls);

const std::map<std::string, std::pair<int, int>> electron_bond_nup_ndn = {
    {"HB", {0, 0}},        {"HEISENBERG", {0, 0}}, {"EXCHANGE", {0, 0}},
    {"ISING", {0, 0}},     {"SZ", {0, 0}},         {"SX", undefined_qns},
    {"SY", undefined_qns}, {"S+", {1, -1}},        {"S-", {-1, 1}},
    {"Cdagup", {1, 0}},    {"Cup", {-1, 0}},       {"Cdagdn", {0, 1}},
    {"Cdn", {0, -1}}};

std::pair<int, int> electron_operator_nup_ndn(BondList const &bonds,
                                              Couplings const &cpls);
std::pair<int, int> tj_operator_nup_ndn(BondList const &bonds,
                                        Couplings const &cpls);

} // namespace hydra::utils
