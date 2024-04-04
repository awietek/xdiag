#ifdef XDIAG_USE_HDF5
#include "file_h5_handler.h"

#include <complex>
#include <vector>

#include <xdiag/io/hdf5/write.h>

namespace xdiag::hdf5 {

using complex = std::complex<double>;

FileH5Handler::FileH5Handler(hid_t file_id, std::string field)
    : file_id_(file_id), field_(field) {}

template <class data_t> void FileH5Handler::operator=(data_t const &data) {
  write(file_id_, field_, data);
}

template void FileH5Handler::operator=(int8_t const &);
template void FileH5Handler::operator=(int16_t const &);
template void FileH5Handler::operator=(int32_t const &);
template void FileH5Handler::operator=(int64_t const &);
template void FileH5Handler::operator=(uint8_t const &);
template void FileH5Handler::operator=(uint16_t const &);
template void FileH5Handler::operator=(uint32_t const &);
template void FileH5Handler::operator=(uint64_t const &);
template void FileH5Handler::operator=(double const &);
template void FileH5Handler::operator=(complex const &);

template void FileH5Handler::operator=(std::vector<int8_t> const &);
template void FileH5Handler::operator=(std::vector<int16_t> const &);
template void FileH5Handler::operator=(std::vector<int32_t> const &);
template void FileH5Handler::operator=(std::vector<int64_t> const &);
template void FileH5Handler::operator=(std::vector<uint8_t> const &);
template void FileH5Handler::operator=(std::vector<uint16_t> const &);
template void FileH5Handler::operator=(std::vector<uint32_t> const &);
template void FileH5Handler::operator=(std::vector<uint64_t> const &);
template void FileH5Handler::operator=(std::vector<double> const &);
template void
FileH5Handler::operator=(std::vector<std::complex<double>> const &);

template void FileH5Handler::operator=(arma::ivec const &);
template void FileH5Handler::operator=(arma::uvec const &);
template void FileH5Handler::operator=(arma::vec const &);
template void FileH5Handler::operator=(arma::cx_vec const &);
template void FileH5Handler::operator=(arma::imat const &);
template void FileH5Handler::operator=(arma::umat const &);
template void FileH5Handler::operator=(arma::mat const &);
template void FileH5Handler::operator=(arma::cx_mat const &);

} // namespace xdiag::hdf5
#endif
