#pragma once
#ifdef XDIAG_USE_HDF5

#include <string>
#include <vector>

#include <xdiag/extern/armadillo/armadillo>
#include <hdf5.h>

namespace xdiag::hdf5 {

template <typename data_t>
void write_scalar(hid_t file_id, std::string field, data_t data);

template <typename data_t>
void write_std_vector(hid_t file_id, std::string field,
                      std::vector<data_t> const &data);

template <typename data_t>
void write_arma_vector(hid_t file_id, std::string field,
                       arma::Col<data_t> const &data);

template <typename data_t>
void write_arma_matrix(hid_t file_id, std::string field,
                       arma::Col<data_t> const &data);

template <typename data_t>
void write(hid_t file_id, std::string field, data_t const &data);

} // namespace xdiag::hdf5
#endif
