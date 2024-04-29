#ifdef XDIAG_USE_HDF5
#include "file_h5_subview.hpp"

#include <complex>
#include <vector>

#include <xdiag/io/hdf5/write.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::hdf5 {

using complex = std::complex<double>;


FileH5Submat::FileH5Submat(hid_t file_id, std::string field, int row_number, int col_number)
    : file_id_(file_id), field_(field), col_number_(col_number), row_number_(row_number) {
    }
template <class data_t> void FileH5Submat::operator=(data_t const &data) {
  if(row_number_<0){
    write_arma_col(file_id_, field_, col_number_, data);
  } else {
    Log.err("other slicing not implemented yet");
  }
}
template void FileH5Submat::operator=(arma::ivec const &);
template void FileH5Submat::operator=(arma::uvec const &);
template void FileH5Submat::operator=(arma::vec const &);
template void FileH5Submat::operator=(arma::cx_vec const &);



FileH5Subcube::FileH5Subcube(hid_t file_id, std::string field, int row_number, int col_number, int slice_number)
    : file_id_(file_id), field_(field), col_number_(col_number), row_number_(row_number), slice_number_(slice_number) {
    }
template <class data_t> void FileH5Subcube::operator=(data_t const &data) {
  if(row_number_<0 && col_number_<0){
    write_arma_slice(file_id_, field_, slice_number_, data);
  } else {
    Log.err("other slicing not implemented yet");
  }
}
template void FileH5Subcube::operator=(arma::imat const &);
template void FileH5Subcube::operator=(arma::umat const &);
template void FileH5Subcube::operator=(arma::mat const &);
template void FileH5Subcube::operator=(arma::cx_mat const &);


} // namespace xdiag::hdf5
#endif
