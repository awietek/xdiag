#include <map>
#include <string>
#include <vector>

#include <xdiag/common.hpp>

namespace xdiag {

inline const std::vector<std::string> known_types = {
    "SdotS",           "Exchange", "SzSz",     "Sz",    "S+",     "S-",
    "ScalarChirality", "Hop",      "Hopup",    "Hopdn", "Cdagup", "Cup",
    "Cdagdn",          "Cdn",      "HubbardU", "Ntot",  "Nup",    "Ndn",
    "tJSzSz",          "tJSdotS",  "Matrix"};

inline const std::vector<std::string> real_types = {
    "SdotS", "Exchange", "SzSz",   "Sz",     "S+",     "S-",  "Hop",
    "Hopup", "Hopdn",    "Cdagup", "Cup",    "Cdagdn", "Cdn", "HubbardU",
    "Ntot",  "Nup",      "Ndn",    "tJSzSz", "tJSdotS"};
inline const std::vector<std::string> cplx_types = {"ScalarChirality"};

inline const std::map<std::string, int64_t> _n_sites_of_type = {
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
    {"tJSzSz", 2},
    {"tJSdotS", 2},
    {"Matrix", undefined}};

bool is_known_type(std::string type);
bool is_real_type(std::string type);
bool is_cplx_type(std::string type);
int64_t n_sites_of_type(std::string type);

std::string known_types_string();

} // namespace xdiag
