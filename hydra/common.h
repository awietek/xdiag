#pragma once

#include <complex>
#include <lila/all.h>

namespace hydra {

inline lila::Logger HydraLog;

using int16 = short;
using int32 = int;
using int64 = long;
using uint16 = unsigned short;
using uint32 = unsigned int;
using uint64 = unsigned long;

using std_bit_t = uint64;
using number_t = int32;

using idx_t = int64;

using scomplex = std::complex<float>;
using complex = std::complex<double>;

} // namespace hydra
