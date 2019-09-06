#ifndef HYDRA_TYPEDEFS_H_
#define HYDRA_TYPEDEFS_H_

#include <complex>

namespace hydra {

  using int16 = short;
  using int32 = int;
  using int64 = long;
  using uint16 = unsigned short;
  using uint32 = unsigned int;
  using uint64 = unsigned long;

  using scomplex = std::complex<float>;
  using complex = std::complex<double>;

}

#endif
