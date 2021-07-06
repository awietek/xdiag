#include "electron_matrix.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/electron/terms/electron_hopping.h>
#include <hydra/models/electron/terms/electron_u.h>

namespace hydra {

template <class bit_t>
lila::Matrix<double>
MatrixReal(BondList const &bonds, Couplings const &couplings,
           Electron<bit_t> const &block_in, Electron<bit_t> const &block_out) {
  assert(block_in == block_out); // only temporary
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  auto mat = lila::Zeros<double>(dim_out, dim_in);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, double val) {
    mat(idx_out, idx_in) += val;
  };

  electron::do_U(couplings, block_in, fill);
  electron::do_hopping<bit_t, double>(bonds, couplings, block_in, fill);

  return mat;
}

template <class bit_t>
lila::Matrix<complex>
MatrixCplx(BondList const &bonds, Couplings const &couplings,
           Electron<bit_t> const &block_in, Electron<bit_t> const &block_out) {
  assert(block_in == block_out);
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  auto mat = lila::Zeros<complex>(dim_out, dim_in);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, complex val) {
    mat(idx_out, idx_in) += val;
  };

  electron::do_U(couplings, block_in, fill);
  electron::do_hopping<bit_t, complex>(bonds, couplings, block_in, fill);

  return mat;
}

template lila::Matrix<double>
MatrixReal<uint16>(BondList const &bonds, Couplings const &couplings,
                   Electron<uint16> const &block_in,
                   Electron<uint16> const &block_out);
template lila::Matrix<double>
MatrixReal<uint32>(BondList const &bonds, Couplings const &couplings,
                   Electron<uint32> const &block_in,
                   Electron<uint32> const &block_out);
template lila::Matrix<double>
MatrixReal<uint64>(BondList const &bonds, Couplings const &couplings,
                   Electron<uint64> const &block_in,
                   Electron<uint64> const &block_out);

template lila::Matrix<complex>
MatrixCplx<uint16>(BondList const &bonds, Couplings const &couplings,
                   Electron<uint16> const &block_in,
                   Electron<uint16> const &block_out);
template lila::Matrix<complex>
MatrixCplx<uint32>(BondList const &bonds, Couplings const &couplings,
                   Electron<uint32> const &block_in,
                   Electron<uint32> const &block_out);
template lila::Matrix<complex>
MatrixCplx<uint64>(BondList const &bonds, Couplings const &couplings,
                   Electron<uint64> const &block_in,
                   Electron<uint64> const &block_out);

} // namespace hydra
