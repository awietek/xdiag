
#include <string>
#include <vector>

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

bool is_known_type(std::string type);
bool is_real_type(std::string type);
bool is_cplx_type(std::string type);

std::string known_types_string();

} // namespace xdiag
