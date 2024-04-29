#pragma once
#ifdef XDIAG_USE_HDF5

#include <string>
#include <vector>

#include <hdf5.h>

#include <xdiag/extern/armadillo/armadillo>

namespace xdiag::hdf5 {

template <typename data_t>
void write_scalar(hid_t file_id, std::string field, data_t data);

template <typename data_t>
void write_std_vector(hid_t file_id, std::string field,
                      std::vector<data_t> const &data);

template <typename data_t>
void write_arma_vector(hid_t file_id, std::string field,
                       arma::Mat<data_t> const &data);

template <typename data_t>
void write_arma_matrix(hid_t file_id, std::string field,
                       arma::Mat<data_t> const &data);

template <typename data_t>
void write_arma_cube(hid_t file_id, std::string field,
                       arma::Cube<data_t> const &data);

template <typename data_t>
void write(hid_t file_id, std::string field, data_t const &data);

// subview operations
template <typename data_t>
void write_arma_col(hid_t file_id, std::string field, int col_number, arma::Col<data_t> const &data);

template <typename data_t>
void write_arma_slice(hid_t file_id, std::string field, int slice_number, arma::Mat<data_t> const &data);


} // namespace xdiag::hdf5
#endif
