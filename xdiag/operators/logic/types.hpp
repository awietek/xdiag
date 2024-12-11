
#include <string>
#include <vector>

namespace xdiag {

inline const std::vector<std::string> known_types = {
    "SDOTS",  "EXCHANGE",        "ISING",    "SZ",      "S+",
    "S-",     "SCALARCHIRALITY", "HOP",      "HOPUP",   "CDAGUP",
    "CUP",    "HOPDN",           "CDAGDN",   "CDN",     "HUBBARDU",
    "NUMBER", "NUMBERUP",        "NUMBERDN", "TJISING", "MATRIX"};

inline const std::vector<std::string> real_types = {
    "SDOTS", "EXCHANGE", "ISING",  "SZ",       "S+",       "S-",
    "HOP",   "HOPUP",    "CDAGUP", "CUP",      "HOPDN",    "CDAGDN",
    "CDN",   "HUBBARDU", "NUMBER", "NUMBERUP", "NUMBERDN", "TJISING"};
inline const std::vector<std::string> cplx_types = {"SCALARCHIRALITY"};

bool is_known_type(std::string type);
bool is_real_type(std::string type);
bool is_cplx_type(std::string type);

} // namespace xdiag