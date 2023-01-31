#include "types.h"

#include <complex>
#include <cstdint>
#include <extern/armadillo/armadillo>

namespace hydra::hdf5 {

using complex = std::complex<double>;

template <> hid_t hdf5_datatype<int8_t>() { return H5T_NATIVE_SCHAR; }
template <> hid_t hdf5_datatype<int16_t>() { return H5T_NATIVE_SHORT; }
template <> hid_t hdf5_datatype<int32_t>() { return H5T_NATIVE_INT; }
template <> hid_t hdf5_datatype<int64_t>() { return H5T_NATIVE_LONG; }
template <> hid_t hdf5_datatype<arma::sword>() { return H5T_NATIVE_LLONG; }

template <> hid_t hdf5_datatype<uint8_t>() { return H5T_NATIVE_UCHAR; }
template <> hid_t hdf5_datatype<uint16_t>() { return H5T_NATIVE_USHORT; }
template <> hid_t hdf5_datatype<uint32_t>() { return H5T_NATIVE_UINT; }
template <> hid_t hdf5_datatype<uint64_t>() { return H5T_NATIVE_ULONG; }
template <> hid_t hdf5_datatype<arma::uword>() { return H5T_NATIVE_ULLONG; }

template <> hid_t hdf5_datatype<double>() { return H5T_NATIVE_DOUBLE; }
template <> hid_t hdf5_datatype<complex>() {
  hid_t memtype = H5Tcreate(H5T_COMPOUND, 2 * sizeof(double));
  H5Tinsert(memtype, "r", 0 * sizeof(double), H5T_NATIVE_DOUBLE);
  H5Tinsert(memtype, "i", 1 * sizeof(double), H5T_NATIVE_DOUBLE);
  return memtype;
}

template <class T> bool hdf5_datatype_mutable() { return false; }
template <> bool hdf5_datatype_mutable<complex>() { return true; }

template bool hdf5_datatype_mutable<int8_t>();
template bool hdf5_datatype_mutable<int16_t>();
template bool hdf5_datatype_mutable<int32_t>();
template bool hdf5_datatype_mutable<int64_t>();
template bool hdf5_datatype_mutable<arma::sword>();

template bool hdf5_datatype_mutable<uint8_t>();
template bool hdf5_datatype_mutable<uint16_t>();
template bool hdf5_datatype_mutable<uint32_t>();
template bool hdf5_datatype_mutable<uint64_t>();
template bool hdf5_datatype_mutable<arma::uword>();

template bool hdf5_datatype_mutable<double>();

} // namespace hydra::hdf5
