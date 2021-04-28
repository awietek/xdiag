#include "electron_matrix.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/electron/terms/electron_hopping.h>
#include <hydra/models/electron/terms/electron_u.h>

namespace hydra {

template <class bit_t>
lila::Matrix<double>
matrix_real(BondList const &bonds, Couplings const &couplings,
            Electron<bit_t> const &block_in, Electron<bit_t> const &block_out) {
  assert(block_in == block_out);  // only temporary
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  auto mat = lila::Zeros<double>(dim_out, dim_in);

  // Hubbard U
  electron::do_U(couplings, block_in,
                 [&mat](idx_t idx, double val) {
                   mat(idx, idx) += val;
                 });

  // electron hopping
  electron::do_hopping<bit_t, double>(
      bonds, couplings, block_in,
      [&mat](idx_t idx_out, idx_t idx_in, double val) {
        mat(idx_out, idx_in) += val;
      });

  return mat;
}
  

template <class bit_t>
lila::Matrix<complex>
matrix_cplx(BondList const &bonds, Couplings const &couplings,
            Electron<bit_t> const &block_in, Electron<bit_t> const &block_out) {
  assert(block_in == block_out);
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  auto mat = lila::Zeros<complex>(dim_out, dim_in);

  // Hubbard U
  electron::do_U(couplings, block_in,
                 [&mat](idx_t idx, double val) {
                   mat(idx, idx) += val;
                 });

  // electron hopping
  electron::do_hopping<bit_t, complex>(
      bonds, couplings, block_in,
      [&mat](idx_t idx_out, idx_t idx_in, complex val) {
        mat(idx_out, idx_in) += val;
      });

  return mat;
}


template lila::Matrix<double>
matrix_real<uint16>(BondList const &bonds, Couplings const &couplings,
                    Electron<uint16> const &block_in,
                    Electron<uint16> const &block_out);
template lila::Matrix<double>
matrix_real<uint32>(BondList const &bonds, Couplings const &couplings,
                    Electron<uint32> const &block_in,
                    Electron<uint32> const &block_out);
template lila::Matrix<double>
matrix_real<uint64>(BondList const &bonds, Couplings const &couplings,
                    Electron<uint64> const &block_in,
                    Electron<uint64> const &block_out);


template lila::Matrix<complex>
matrix_cplx<uint16>(BondList const &bonds, Couplings const &couplings,
                    Electron<uint16> const &block_in,
                    Electron<uint16> const &block_out);
template lila::Matrix<complex>
matrix_cplx<uint32>(BondList const &bonds, Couplings const &couplings,
                    Electron<uint32> const &block_in,
                    Electron<uint32> const &block_out);
template lila::Matrix<complex>
matrix_cplx<uint64>(BondList const &bonds, Couplings const &couplings,
                    Electron<uint64> const &block_in,
                    Electron<uint64> const &block_out);

} // namespace hydra
