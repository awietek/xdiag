#pragma once
#ifdef HYDRA_USE_HDF5

#include <string>
#include <vector>

#include <hdf5.h>
#include <extern/armadillo/armadillo>

namespace hydra::hdf5 {

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

} // namespace hydra::hdf5
#endif
