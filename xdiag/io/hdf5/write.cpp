// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#ifdef XDIAG_USE_HDF5
#include "write.hpp"

#include <complex>
#include <cstdint>
#include <vector>

#include <xdiag/io/hdf5/types.hpp>
#include <xdiag/io/hdf5/utils.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::hdf5 {

using complex = std::complex<double>;

template <typename data_t>
void write_scalar(hid_t file_id, std::string field, data_t data) try {
  std::vector<hid_t> groups = create_groups(file_id, field);
  hid_t group = (groups.size() == 0) ? file_id : groups[groups.size() - 1];
  if (group == H5I_INVALID_HID) {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error creating groups for field \"{}\"", field));
  }

  hid_t datatype = hdf5_datatype<data_t>();
  hid_t dataspace = H5Screate(H5S_SCALAR);
  std::string name = dataset_name(field);

  hid_t dataset;
  if (H5Lexists(file_id, field.c_str(), H5P_DEFAULT)) {
    dataset = H5Dopen(file_id, field.c_str(), H5P_DEFAULT);
  } else {
    dataset = H5Dcreate(group, name.c_str(), datatype, dataspace, H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT);
  }

  if (dataset != H5I_INVALID_HID) {
    hid_t status =
        H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
    if (status < 0) {
      XDIAG_THROW(fmt::format(
          "Error in xdiag hdf5: could not write data for field \"{}\"", field));
    }
  } else {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error creating dataset for field \"{}\"", field));
  }

  H5Dclose(dataset);
  H5Sclose(dataspace);
  if (hdf5_datatype_mutable<data_t>()) {
    H5Tclose(datatype);
  }
  close_groups(groups);
}
XDIAG_CATCH

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
                      std::vector<data_t> const &data) try {

  std::vector<hid_t> groups = create_groups(file_id, field);
  hid_t group = (groups.size() == 0) ? file_id : groups[groups.size() - 1];
  if (group == H5I_INVALID_HID) {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error creating groups for field \"{}\"", field));
  }
  hid_t datatype = hdf5_datatype<data_t>();
  hsize_t dims[1];
  dims[0] = data.size();
  hid_t dataspace = H5Screate_simple(1, dims, nullptr);
  std::string name = dataset_name(field);

  hid_t dataset;
  if (H5Lexists(file_id, field.c_str(), H5P_DEFAULT)) {
    dataset = H5Dopen(file_id, field.c_str(), H5P_DEFAULT);
  } else {
    dataset = H5Dcreate(group, name.c_str(), datatype, dataspace, H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT);
  }

  if (dataset != H5I_INVALID_HID) {
    hid_t status =
        H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    if (status < 0) {
      XDIAG_THROW(fmt::format(
          "Error in xdiag hdf5: could not write data for field \"{}\"", field));
    }
  } else {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error creating dataset for field \"{}\"", field));
  }
  H5Dclose(dataset);
  H5Sclose(dataspace);
  if (hdf5_datatype_mutable<data_t>()) {
    H5Tclose(datatype);
  }
  close_groups(groups);
}
XDIAG_CATCH

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
                       arma::Col<data_t> const &data) try {

  std::vector<hid_t> groups = create_groups(file_id, field);
  hid_t group = (groups.size() == 0) ? file_id : groups[groups.size() - 1];
  if (group == H5I_INVALID_HID) {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error creating groups for field \"{}\"", field));
  }
  hid_t datatype = hdf5_datatype<data_t>();
  hsize_t dims[1];
  dims[0] = data.n_elem;
  hid_t dataspace = H5Screate_simple(1, dims, nullptr);
  std::string name = dataset_name(field);

  hid_t dataset;
  if (H5Lexists(file_id, field.c_str(), H5P_DEFAULT)) {
    dataset = H5Dopen(file_id, field.c_str(), H5P_DEFAULT);
  } else {
    dataset = H5Dcreate(group, name.c_str(), datatype, dataspace, H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT);
  }

  if (dataset != H5I_INVALID_HID) {
    hid_t status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            data.memptr());
    if (status < 0) {
      XDIAG_THROW(fmt::format(
          "Error in xdiag hdf5: could not write data for field \"{}\"", field));
    }
  } else {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error creating dataset for field \"{}\"", field));
  }
  H5Dclose(dataset);
  H5Sclose(dataspace);
  if (hdf5_datatype_mutable<data_t>()) {
    H5Tclose(datatype);
  }
  close_groups(groups);
}
XDIAG_CATCH

template void write_arma_vector(hid_t, std::string,
                                arma::Col<arma::sword> const &);
template void write_arma_vector(hid_t, std::string,
                                arma::Col<arma::uword> const &);
template void write_arma_vector(hid_t, std::string, arma::Col<double> const &);
template void write_arma_vector(hid_t, std::string, arma::Col<complex> const &);

template <typename data_t>
void write_arma_matrix(hid_t file_id, std::string field,
                       arma::Mat<data_t> const &data) try {

  std::vector<hid_t> groups = create_groups(file_id, field);
  hid_t group = (groups.size() == 0) ? file_id : groups[groups.size() - 1];
  if (group == H5I_INVALID_HID) {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error creating groups for field \"{}\"", field));
  }
  hid_t datatype = hdf5_datatype<data_t>();
  hsize_t dims[2];
  dims[1] = data.n_rows;
  dims[0] = data.n_cols;

  hid_t dataspace = H5Screate_simple(2, dims, nullptr);
  std::string name = dataset_name(field);

  hid_t dataset;
  if (H5Lexists(file_id, field.c_str(), H5P_DEFAULT)) {
    dataset = H5Dopen(file_id, field.c_str(), H5P_DEFAULT);
  } else {
    dataset = H5Dcreate(group, name.c_str(), datatype, dataspace, H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT);
  }

  if (dataset != H5I_INVALID_HID) {

    hid_t status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            data.memptr());

    if (status < 0) {
      XDIAG_THROW(fmt::format(
          "Error in xdiag hdf5: could not write data for field \"{}\"", field));
    }
  } else {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error creating dataset for field \"{}\"", field));
  }
  H5Dclose(dataset);
  H5Sclose(dataspace);
  if (hdf5_datatype_mutable<data_t>()) {
    H5Tclose(datatype);
  }
  close_groups(groups);
}
XDIAG_CATCH

template void write_arma_matrix(hid_t, std::string,
                                arma::Mat<arma::sword> const &);
template void write_arma_matrix(hid_t, std::string,
                                arma::Mat<arma::uword> const &);
template void write_arma_matrix(hid_t, std::string, arma::Mat<double> const &);
template void write_arma_matrix(hid_t, std::string, arma::Mat<complex> const &);

// arma cubes ==================================================================

template <typename data_t>
void write_arma_cube(hid_t file_id, std::string field,
                     arma::Cube<data_t> const &data) try {

  std::vector<hid_t> groups = create_groups(file_id, field);
  hid_t group = (groups.size() == 0) ? file_id : groups[groups.size() - 1];
  if (group == H5I_INVALID_HID) {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error creating groups for field \"{}\"", field));
  }
  hid_t datatype = hdf5_datatype<data_t>();
  hsize_t dims[3];
  dims[2] = data.n_rows;
  dims[1] = data.n_cols;
  dims[0] = data.n_slices;

  hid_t dataspace = H5Screate_simple(3, dims, nullptr);
  std::string name = dataset_name(field);

  hid_t dataset;
  if (H5Lexists(file_id, field.c_str(), H5P_DEFAULT)) {
    dataset = H5Dopen(file_id, field.c_str(), H5P_DEFAULT);
  } else {
    dataset = H5Dcreate(group, name.c_str(), datatype, dataspace, H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT);
  }

  if (dataset != H5I_INVALID_HID) {

    hid_t status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            data.memptr());

    if (status < 0) {
      XDIAG_THROW(fmt::format(
          "Error in xdiag hdf5: could not write data for field \"{}\"", field));
    }
  } else {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error creating dataset for field \"{}\"", field));
  }
  H5Dclose(dataset);
  H5Sclose(dataspace);
  if (hdf5_datatype_mutable<data_t>()) {
    H5Tclose(datatype);
  }
  close_groups(groups);
}
XDIAG_CATCH

template void write_arma_cube(hid_t, std::string,
                              arma::Cube<arma::sword> const &);
template void write_arma_cube(hid_t, std::string,
                              arma::Cube<arma::uword> const &);
template void write_arma_cube(hid_t, std::string, arma::Cube<double> const &);
template void write_arma_cube(hid_t, std::string, arma::Cube<complex> const &);

// submatrix operations ========================================================
template <typename data_t>
void write_arma_col(hid_t file_id, std::string field, int col_number,
                    arma::Col<data_t> const &data) try {

  std::vector<hid_t> groups = create_groups(file_id, field);
  hid_t group = (groups.size() == 0) ? file_id : groups[groups.size() - 1];
  if (group == H5I_INVALID_HID) {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error creating groups for field \"{}\"", field));
  }
  if (!H5Lexists(file_id, field.c_str(), H5P_DEFAULT)) {
    XDIAG_THROW(fmt::format("data set must exist to append col."));
  }
  // load the old dataset
  hid_t dataset = H5Dopen(file_id, field.c_str(), H5P_DEFAULT);

  // reshape the col into an arma mat
  arma::Mat<data_t> data_r(data.n_elem, 1);
  data_r.col(0) = data;

  // create the mem space for the new dataset
  hsize_t rank = 2;
  hsize_t dims[2];
  dims[1] = data_r.n_rows;
  dims[0] = data_r.n_cols;
  hid_t mem_space = H5Screate_simple(rank, dims, NULL); // shape of new data
  hid_t file_space = H5Dget_space(dataset);             // shape of disk data

  // create the hyperslab
  hsize_t start[2];
  start[1] = 0;
  start[0] = col_number;

  hsize_t count[2];
  count[1] = data_r.n_rows;
  count[0] = 1;

  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);

  hid_t datatype = hdf5_datatype<data_t>();
  std::string name = dataset_name(field);

  if (dataset != H5I_INVALID_HID) {
    hid_t status = H5Dwrite(dataset, datatype, mem_space, file_space,
                            H5P_DEFAULT, data_r.memptr());

    if (status < 0) {
      XDIAG_THROW(fmt::format(
          "Error in xdiag hdf5: could not append col data for field \"{}\"",
          field));
    }
  } else {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error opening dataset for field \"{}\"", field));
  }
  H5Sclose(file_space);
  H5Sclose(mem_space);
  H5Dclose(dataset);

  if (hdf5_datatype_mutable<data_t>()) {
    H5Tclose(datatype);
  }
  close_groups(groups);
}
XDIAG_CATCH

template void write_arma_col(hid_t, std::string, int col_number,
                             arma::Col<arma::sword> const &);
template void write_arma_col(hid_t, std::string, int col_number,
                             arma::Col<arma::uword> const &);
template void write_arma_col(hid_t, std::string, int col_number,
                             arma::Col<double> const &);
template void write_arma_col(hid_t, std::string, int col_number,
                             arma::Col<complex> const &);

// updata and dynamically allocat the memspace?
template <typename data_t>
void write_arma_slice(hid_t file_id, std::string field, int slice_number,
                      arma::Mat<data_t> const &data) try {

  std::vector<hid_t> groups = create_groups(file_id, field);
  hid_t group = (groups.size() == 0) ? file_id : groups[groups.size() - 1];
  if (group == H5I_INVALID_HID) {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error creating groups for field \"{}\"", field));
  }
  if (!H5Lexists(file_id, field.c_str(), H5P_DEFAULT)) {
    XDIAG_THROW(fmt::format("data set must exist to append col."));
  }
  // load the old dataset
  hid_t dataset = H5Dopen(file_id, field.c_str(), H5P_DEFAULT);

  // reshape the col into an arma mat
  arma::Cube<data_t> data_r(data.n_rows, data.n_cols, 1);
  data_r.slice(0) = data;

  // create the mem space for the new dataset
  hsize_t rank = 3;
  hsize_t dims[3];
  dims[2] = data_r.n_rows;
  dims[1] = data_r.n_cols;
  dims[0] = data_r.n_slices;
  hid_t mem_space = H5Screate_simple(rank, dims, NULL); // shape of new data
  hid_t file_space = H5Dget_space(dataset);             // shape of disk data

  // create the hyperslab
  hsize_t start[3];
  start[2] = 0;
  start[1] = 0;
  start[0] = slice_number;

  hsize_t count[3];
  count[2] = data_r.n_rows;
  count[1] = data_r.n_cols;
  count[0] = 1;

  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);

  hid_t datatype = hdf5_datatype<data_t>();
  std::string name = dataset_name(field);

  if (dataset != H5I_INVALID_HID) {
    hid_t status = H5Dwrite(dataset, datatype, mem_space, file_space,
                            H5P_DEFAULT, data_r.memptr());

    if (status < 0) {
      XDIAG_THROW(fmt::format(
          "Error in xdiag hdf5: could not append col data for field \"{}\"",
          field));
    }
  } else {
    XDIAG_THROW(fmt::format(
        "Error in xdiag hdf5: error opening dataset for field \"{}\"", field));
  }
  H5Sclose(file_space);
  H5Sclose(mem_space);
  H5Dclose(dataset);

  if (hdf5_datatype_mutable<data_t>()) {
    H5Tclose(datatype);
  }
  close_groups(groups);
}
XDIAG_CATCH

template void write_arma_slice(hid_t, std::string, int slice_number,
                               arma::Mat<arma::sword> const &);
template void write_arma_slice(hid_t, std::string, int slice_number,
                               arma::Mat<arma::uword> const &);
template void write_arma_slice(hid_t, std::string, int slice_number,
                               arma::Mat<double> const &);
template void write_arma_slice(hid_t, std::string, int slice_number,
                               arma::Mat<complex> const &);

template <>
void write(hid_t file_id, std::string field, int8_t const &data) try {
  write_scalar(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field, int16_t const &data) try {
  write_scalar(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field, int32_t const &data) try {
  write_scalar(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field, int64_t const &data) try {
  write_scalar(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field, uint8_t const &data) try {
  write_scalar(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field, uint16_t const &data) try {
  write_scalar(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field, uint32_t const &data) try {
  write_scalar(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field, uint64_t const &data) try {
  write_scalar(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field, double const &data) try {
  write_scalar(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field, complex const &data) try {
  write_scalar(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           std::vector<int8_t> const &data) try {
  write_std_vector(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           std::vector<int16_t> const &data) try {
  write_std_vector(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           std::vector<int32_t> const &data) try {
  write_std_vector(file_id, field, data);
}
XDIAG_CATCH
template <>
void write(hid_t file_id, std::string field,
           std::vector<int64_t> const &data) try {
  write_std_vector(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           std::vector<uint8_t> const &data) try {
  write_std_vector(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           std::vector<uint16_t> const &data) try {
  write_std_vector(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           std::vector<uint32_t> const &data) try {
  write_std_vector(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           std::vector<uint64_t> const &data) try {
  write_std_vector(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           std::vector<double> const &data) try {
  write_std_vector(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           std::vector<complex> const &data) try {
  write_std_vector(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           arma::Col<arma::sword> const &data) try {
  write_arma_vector(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           arma::Col<arma::uword> const &data) try {
  write_arma_vector(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           arma::Col<double> const &data) try {
  write_arma_vector(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           arma::Col<complex> const &data) try {
  write_arma_vector(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           arma::Mat<arma::sword> const &data) try {
  write_arma_matrix(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           arma::Mat<arma::uword> const &data) try {
  write_arma_matrix(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           arma::Mat<double> const &data) try {
  write_arma_matrix(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           arma::Mat<complex> const &data) try {
  write_arma_matrix(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           arma::Cube<arma::sword> const &data) try {
  write_arma_cube(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           arma::Cube<arma::uword> const &data) try {
  write_arma_cube(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           arma::Cube<double> const &data) try {
  write_arma_cube(file_id, field, data);
}
XDIAG_CATCH

template <>
void write(hid_t file_id, std::string field,
           arma::Cube<complex> const &data) try {
  write_arma_cube(file_id, field, data);
}
XDIAG_CATCH

} // namespace xdiag::hdf5
#endif
