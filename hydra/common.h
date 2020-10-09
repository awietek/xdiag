#ifndef HYDRA_COMMON_H_
#define HYDRA_COMMON_H_

#include <complex>

namespace hydra {

  using int16 = short;
  using int32 = int;
  using int64 = long;
  using uint16 = unsigned short;
  using uint32 = unsigned int;
  using uint64 = unsigned long;
  
  using std_bit_t = uint64;
  using number_t = int32;

  using std_idx_t = uint64;
  
  using scomplex = std::complex<float>;
  using complex = std::complex<double>;
  
}

#endif
