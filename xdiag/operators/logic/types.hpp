// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <map>
#include <string>
#include <vector>

#include <xdiag/common.hpp>

namespace xdiag {

inline const std::vector<std::string> known_types = {
    "Id",    "SdotS",           "Exchange", "SzSz",       "Sz",     "S+",
    "S-",    "ScalarChirality", "Hop",      "Hopup",      "Hopdn",  "Cdagup",
    "Cup",   "Cdagdn",          "Cdn",      "HubbardU",   "Ntot",   "Nup",
    "Ndn",   "Nupdn",           "NtotNtot", "NupdnNupdn", "tJSzSz", "tJSdotS",
    "Matrix", "Nhup", "Nhdn", "NupNdn", "NhupNdn", "NupNhdn", "NhupNhdn",
    "Nhtot", "NhtotNhtot", "NhtotNtot", "NupNup", "NdnNdn", "NhupNup", "NhdnNdn",
    "NhupNhup", "NhdnNhdn"};

inline const std::vector<std::string> real_types = {
    "Id",     "SdotS",    "Exchange",   "SzSz",   "Sz",     "S+",
    "S-",     "Hop",      "Hopup",      "Hopdn",  "Cdagup", "Cup",
    "Cdagdn", "Cdn",      "HubbardU",   "Ntot",   "Nup",    "Ndn",
    "Nupdn",  "NtotNtot", "NupdnNupdn", "tJSzSz", "tJSdotS", "Nhup", "Nhdn",
    "NupNdn", "NhupNdn", "NupNhdn", "NhupNhdn", "Nhtot", "NhtotNhtot", "NhtotNtot", "NupNup",
    "NdnNdn", "NhupNup", "NhdnNdn", "NhupNhup", "NhdnNhdn"};
inline const std::vector<std::string> cplx_types = {"ScalarChirality"};

inline const std::map<std::string, int64_t> _nsites_of_type = {
    {"Id", undefined},
    {"SdotS", 2},
    {"Exchange", 2},
    {"SzSz", 2},
    {"Sz", 1},
    {"S+", 1},
    {"S-", 1},
    {"ScalarChirality", 3},
    {"Hop", 2},
    {"Hopup", 2},
    {"Hopdn", 2},
    {"Cdagup", 1},
    {"Cup", 1},
    {"Cdagdn", 1},
    {"Cdn", 1},
    {"HubbardU", undefined},
    {"Ntot", 1},
    {"Nup", 1},
    {"Ndn", 1},
    {"Nupdn", 1},
    {"NtotNtot", 2},
    {"NupdnNupdn", 2},
    {"tJSzSz", 2},
    {"tJSdotS", 2},
    {"Matrix", undefined},
    {"Nhup", 1},
    {"Nhdn", 1},
    {"NupNdn", 2},
    {"NhupNdn", 2},
    {"NupNhdn", 2},
    {"NhupNhdn", 2},
    {"Nhtot", 1},
    {"NhtotNhtot", 2},
    {"NhtotNtot", 2},
    {"NupNup", 2},
    {"NdnNdn", 2},
    {"NhupNup", 2},
    {"NhdnNdn", 2},
    {"NhupNhup", 2},
    {"NhdnNhdn", 2}
};

bool is_known_type(std::string type);
bool is_real_type(std::string type);
bool is_cplx_type(std::string type);
int64_t nsites_of_type(std::string type);

std::string known_types_string();

} // namespace xdiag
