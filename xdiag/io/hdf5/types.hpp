#pragma once
#ifdef XDIAG_USE_HDF5

#include <string>
#include <vector>

#include <hdf5.h>

namespace xdiag::hdf5 {

template <class T> hid_t hdf5_datatype();
template <class T> bool hdf5_datatype_mutable();

} // namespace xdiag::hdf5
#endif
