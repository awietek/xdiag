#ifdef XDIAG_USE_HDF5
#include "write.h"

#include <complex>
#include <cstdint>
#include <vector>

#include <xdiag/io/hdf5/types.h>
#include <xdiag/io/hdf5/utils.h>
#include <xdiag/utils/logger.h>

namespace xdiag::hdf5 {

using complex = std::complex<double>;

template <typename data_t>
void write_scalar(hid_t file_id, std::string field, data_t data) {
  std::vector<hid_t> groups = create_groups(file_id, field);
  hid_t group = (groups.size() == 0) ? file_id : groups[groups.size() - 1];
  if (group == H5I_INVALID_HID) {
    Log.err("Error in xdiag hdf5: error creating groups for field \"{}\"",
            field);
  }

  hid_t datatype = hdf5_datatype<data_t>();
  hid_t dataspace = H5Screate(H5S_SCALAR);
  std::string name = dataset_name(field);
  hid_t dataset = H5Dcreate(group, name.c_str(), datatype, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset != H5I_INVALID_HID) {
    hid_t status =
        H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
    if (status < 0) {
      Log.err("Error in xdiag hdf5: could not write data for field \"{}\"",
              field);
    }
  } else {
    Log.err("Error in xdiag hdf5: error creating dataset for field \"{}\"",
            field);
  }

  H5Dclose(dataset);
  H5Sclose(dataspace);
  if (hdf5_datatype_mutable<data_t>()) {
    H5Tclose(datatype);
  }
  close_groups(groups);
}

template void write_scalar(hid_t, std::string, int8_t);
template void write_scalar(hid_t, std::string, int16_t);
template void write_scalar(hid_t, std::string, int32_t);
template void write_scalar(hid_t, std::string, int64_t);
template void write_scalar(hid_t, std::string, uint8_t);
template void write_scalar(hid_t, std::string, uint16_t);
template void write_scalar(hid_t, std::string, uint32_t);
template void write_scalar(hid_t, std::string, uint64_t);
template void write_scalar(hid_t, std::string, double);
template void write_scalar(hid_t, std::string, complex);

template <typename data_t>
void write_std_vector(hid_t file_id, std::string field,
                      std::vector<data_t> const &data) {

  std::vector<hid_t> groups = create_groups(file_id, field);
  hid_t group = (groups.size() == 0) ? file_id : groups[groups.size() - 1];
  if (group == H5I_INVALID_HID) {
    Log.err("Error in xdiag hdf5: error creating groups for field \"{}\"",
            field);
  }
  hid_t datatype = hdf5_datatype<data_t>();
  hsize_t dims[1];
  dims[0] = data.size();
  hid_t dataspace = H5Screate_simple(1, dims, nullptr);
  std::string name = dataset_name(field);
  hid_t dataset = H5Dcreate(group, name.c_str(), datatype, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset != H5I_INVALID_HID) {
    hid_t status =
        H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    if (status < 0) {
      Log.err("Error in xdiag hdf5: could not write data for field \"{}\"",
              field);
    }
  } else {
    Log.err("Error in xdiag hdf5: error creating dataset for field \"{}\"",
            field);
  }
  H5Dclose(dataset);
  H5Sclose(dataspace);
  if (hdf5_datatype_mutable<data_t>()) {
    H5Tclose(datatype);
  }
  close_groups(groups);
}
template void write_std_vector(hid_t, std::string, std::vector<int8_t> const &);
template void write_std_vector(hid_t, std::string,
                               std::vector<int16_t> const &);
template void write_std_vector(hid_t, std::string,
                               std::vector<int32_t> const &);
template void write_std_vector(hid_t, std::string,
                               std::vector<int64_t> const &);
template void write_std_vector(hid_t, std::string,
                               std::vector<uint8_t> const &);
template void write_std_vector(hid_t, std::string,
                               std::vector<uint16_t> const &);
template void write_std_vector(hid_t, std::string,
                               std::vector<uint32_t> const &);
template void write_std_vector(hid_t, std::string,
                               std::vector<uint64_t> const &);
template void write_std_vector(hid_t, std::string, std::vector<double> const &);
template void write_std_vector(hid_t, std::string,
                               std::vector<complex> const &);

template <typename data_t>
void write_arma_vector(hid_t file_id, std::string field,
                       arma::Col<data_t> const &data) {

  std::vector<hid_t> groups = create_groups(file_id, field);
  hid_t group = (groups.size() == 0) ? file_id : groups[groups.size() - 1];
  if (group == H5I_INVALID_HID) {
    Log.err("Error in xdiag hdf5: error creating groups for field \"{}\"",
            field);
  }
  hid_t datatype = hdf5_datatype<data_t>();
  hsize_t dims[1];
  dims[0] = data.n_elem;
  hid_t dataspace = H5Screate_simple(1, dims, nullptr);
  std::string name = dataset_name(field);
  hid_t dataset = H5Dcreate(group, name.c_str(), datatype, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset != H5I_INVALID_HID) {
    hid_t status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            data.memptr());
    if (status < 0) {
      Log.err("Error in xdiag hdf5: could not write data for field \"{}\"",
              field);
    }
  } else {
    Log.err("Error in xdiag hdf5: error creating dataset for field \"{}\"",
            field);
  }
  H5Dclose(dataset);
  H5Sclose(dataspace);
  if (hdf5_datatype_mutable<data_t>()) {
    H5Tclose(datatype);
  }
  close_groups(groups);
}

template void write_arma_vector(hid_t, std::string,
                                arma::Col<arma::sword> const &);
template void write_arma_vector(hid_t, std::string,
                                arma::Col<arma::uword> const &);
template void write_arma_vector(hid_t, std::string, arma::Col<double> const &);
template void write_arma_vector(hid_t, std::string, arma::Col<complex> const &);

template <typename data_t>
void write_arma_matrix(hid_t file_id, std::string field,
                       arma::Mat<data_t> const &data) {

  std::vector<hid_t> groups = create_groups(file_id, field);
  hid_t group = (groups.size() == 0) ? file_id : groups[groups.size() - 1];
  if (group == H5I_INVALID_HID) {
    Log.err("Error in xdiag hdf5: error creating groups for field \"{}\"",
            field);
  }
  hid_t datatype = hdf5_datatype<data_t>();
  hsize_t dims[2];
  dims[0] = data.n_rows;
  dims[1] = data.n_cols;
  
  hid_t dataspace = H5Screate_simple(2, dims, nullptr);
  std::string name = dataset_name(field);
  hid_t dataset = H5Dcreate(group, name.c_str(), datatype, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset != H5I_INVALID_HID) {

    arma::Mat<data_t> data_trans = data.st(); // IIIIIEEEEEK
    hid_t status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            data_trans.memptr());

    if (status < 0) {
      Log.err("Error in xdiag hdf5: could not write data for field \"{}\"",
              field);
    }
  } else {
    Log.err("Error in xdiag hdf5: error creating dataset for field \"{}\"",
            field);
  }
  H5Dclose(dataset);
  H5Sclose(dataspace);
  if (hdf5_datatype_mutable<data_t>()) {
    H5Tclose(datatype);
  }
  close_groups(groups);
}

template void write_arma_matrix(hid_t, std::string,
                                arma::Mat<arma::sword> const &);
template void write_arma_matrix(hid_t, std::string,
                                arma::Mat<arma::uword> const &);
template void write_arma_matrix(hid_t, std::string, arma::Mat<double> const &);
template void write_arma_matrix(hid_t, std::string, arma::Mat<complex> const &);

template <> void write(hid_t file_id, std::string field, int8_t const &data) {
  write_scalar(file_id, field, data);
}
template <> void write(hid_t file_id, std::string field, int16_t const &data) {
  write_scalar(file_id, field, data);
}
template <> void write(hid_t file_id, std::string field, int32_t const &data) {
  write_scalar(file_id, field, data);
}
template <> void write(hid_t file_id, std::string field, int64_t const &data) {
  write_scalar(file_id, field, data);
}
template <> void write(hid_t file_id, std::string field, uint8_t const &data) {
  write_scalar(file_id, field, data);
}
template <> void write(hid_t file_id, std::string field, uint16_t const &data) {
  write_scalar(file_id, field, data);
}
template <> void write(hid_t file_id, std::string field, uint32_t const &data) {
  write_scalar(file_id, field, data);
}
template <> void write(hid_t file_id, std::string field, uint64_t const &data) {
  write_scalar(file_id, field, data);
}
template <> void write(hid_t file_id, std::string field, double const &data) {
  write_scalar(file_id, field, data);
}
template <> void write(hid_t file_id, std::string field, complex const &data) {
  write_scalar(file_id, field, data);
}

template <>
void write(hid_t file_id, std::string field, std::vector<int8_t> const &data) {
  write_std_vector(file_id, field, data);
}
template <>
void write(hid_t file_id, std::string field, std::vector<int16_t> const &data) {
  write_std_vector(file_id, field, data);
}
template <>
void write(hid_t file_id, std::string field, std::vector<int32_t> const &data) {
  write_std_vector(file_id, field, data);
}
template <>
void write(hid_t file_id, std::string field, std::vector<int64_t> const &data) {
  write_std_vector(file_id, field, data);
}
template <>
void write(hid_t file_id, std::string field, std::vector<uint8_t> const &data) {
  write_std_vector(file_id, field, data);
}
template <>
void write(hid_t file_id, std::string field,
           std::vector<uint16_t> const &data) {
  write_std_vector(file_id, field, data);
}
template <>
void write(hid_t file_id, std::string field,
           std::vector<uint32_t> const &data) {
  write_std_vector(file_id, field, data);
}
template <>
void write(hid_t file_id, std::string field,
           std::vector<uint64_t> const &data) {
  write_std_vector(file_id, field, data);
}
template <>
void write(hid_t file_id, std::string field, std::vector<double> const &data) {
  write_std_vector(file_id, field, data);
}
template <>
void write(hid_t file_id, std::string field, std::vector<complex> const &data) {
  write_std_vector(file_id, field, data);
}

template <>
void write(hid_t file_id, std::string field,
           arma::Col<arma::sword> const &data) {
  write_arma_vector(file_id, field, data);
}
template <>
void write(hid_t file_id, std::string field,
           arma::Col<arma::uword> const &data) {
  write_arma_vector(file_id, field, data);
}

template <>
void write(hid_t file_id, std::string field, arma::Col<double> const &data) {
  write_arma_vector(file_id, field, data);
}
template <>
void write(hid_t file_id, std::string field, arma::Col<complex> const &data) {
  write_arma_vector(file_id, field, data);
}

template <>
void write(hid_t file_id, std::string field,
           arma::Mat<arma::sword> const &data) {
  write_arma_matrix(file_id, field, data);
}
template <>
void write(hid_t file_id, std::string field,
           arma::Mat<arma::uword> const &data) {
  write_arma_matrix(file_id, field, data);
}

template <>
void write(hid_t file_id, std::string field, arma::Mat<double> const &data) {
  write_arma_matrix(file_id, field, data);
}
template <>
void write(hid_t file_id, std::string field, arma::Mat<complex> const &data) {
  write_arma_matrix(file_id, field, data);
}

} // namespace xdiag::hdf5
#endif
