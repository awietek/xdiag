#pragma once
#ifdef HYDRA_USE_HDF5

#include <string>
#include <vector>

#include <hdf5.h>

namespace hydra::hdf5 {

template <class T> hid_t hdf5_datatype();
template <class T> bool hdf5_datatype_mutable();

} // namespace hydra::hdf5
#endif
