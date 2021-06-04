#include "spinhalf_matrix.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/spinhalf/terms/spinhalf_exchange.h>
#include <hydra/models/spinhalf/terms/spinhalf_ising.h>

namespace hydra {

template <class bit_t>
lila::Matrix<double>
matrix_real(BondList const &bonds, Couplings const &couplings,
            Spinhalf<bit_t> const &block_in, Spinhalf<bit_t> const &block_out) {
  assert(block_in == block_out); // only temporary
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  auto mat = lila::Zeros<double>(dim_out, dim_in);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, double val) {
    mat(idx_out, idx_in) += val;
  };

  spinhalfdetail::do_ising(bonds, couplings, block_in, fill);
  spinhalfdetail::do_exchange<bit_t, complex>(bonds, couplings, block_in, fill);

  return mat;
}

template <class bit_t>
lila::Matrix<complex>
matrix_cplx(BondList const &bonds, Couplings const &couplings,
            Spinhalf<bit_t> const &block_in, Spinhalf<bit_t> const &block_out) {
  assert(block_in == block_out);
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  auto mat = lila::Zeros<complex>(dim_out, dim_in);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, complex val) {
    mat(idx_out, idx_in) += val;
  };

  spinhalfdetail::do_ising(bonds, couplings, block_in, fill);
  spinhalfdetail::do_exchange<bit_t, complex>(bonds, couplings, block_in, fill);

  return mat;
}

template lila::Matrix<double>
matrix_real<uint16>(BondList const &bonds, Couplings const &couplings,
                    Spinhalf<uint16> const &block_in,
                    Spinhalf<uint16> const &block_out);
template lila::Matrix<double>
matrix_real<uint32>(BondList const &bonds, Couplings const &couplings,
                    Spinhalf<uint32> const &block_in,
                    Spinhalf<uint32> const &block_out);
template lila::Matrix<double>
matrix_real<uint64>(BondList const &bonds, Couplings const &couplings,
                    Spinhalf<uint64> const &block_in,
                    Spinhalf<uint64> const &block_out);

template lila::Matrix<complex>
matrix_cplx<uint16>(BondList const &bonds, Couplings const &couplings,
                    Spinhalf<uint16> const &block_in,
                    Spinhalf<uint16> const &block_out);
template lila::Matrix<complex>
matrix_cplx<uint32>(BondList const &bonds, Couplings const &couplings,
                    Spinhalf<uint32> const &block_in,
                    Spinhalf<uint32> const &block_out);
template lila::Matrix<complex>
matrix_cplx<uint64>(BondList const &bonds, Couplings const &couplings,
                    Spinhalf<uint64> const &block_in,
                    Spinhalf<uint64> const &block_out);

} // namespace hydra
