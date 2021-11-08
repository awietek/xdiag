#include "spinhalf_matrix.h"

#include <hydra/blocks/spinhalf/terms/spinhalf_exchange.h>
#include <hydra/blocks/spinhalf/terms/spinhalf_ising.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <class bit_t>
lila::Matrix<double>
MatrixReal(BondList const &bonds, Couplings const &couplings,
           Spinhalf<bit_t> const &block_in, Spinhalf<bit_t> const &block_out) {
  using namespace hydra::terms::spinhalf;

  assert(block_in == block_out); // only temporary

  utils::check_operator_real(bonds, couplings,
                             "construct real Spinhalf matrix");

  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  auto mat = lila::Zeros<double>(dim_out, dim_in);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, double val) {
    mat(idx_out, idx_in) += val;
  };

  do_ising(bonds, couplings, block_in, fill);
  do_exchange<bit_t, double>(bonds, couplings, block_in, fill);

  return mat;
}

template <class bit_t>
lila::Matrix<complex>
MatrixCplx(BondList const &bonds, Couplings const &couplings,
           Spinhalf<bit_t> const &block_in, Spinhalf<bit_t> const &block_out) {
  using namespace terms::spinhalf;

  assert(block_in == block_out);
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  auto mat = lila::Zeros<complex>(dim_out, dim_in);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, complex val) {
    mat(idx_out, idx_in) += val;
  };

  do_ising(bonds, couplings, block_in, fill);
  do_exchange<bit_t, complex>(bonds, couplings, block_in, fill);

  return mat;
}

template lila::Matrix<double>
MatrixReal<uint16_t>(BondList const &bonds, Couplings const &couplings,
                     Spinhalf<uint16_t> const &block_in,
                     Spinhalf<uint16_t> const &block_out);
template lila::Matrix<double>
MatrixReal<uint32_t>(BondList const &bonds, Couplings const &couplings,
                     Spinhalf<uint32_t> const &block_in,
                     Spinhalf<uint32_t> const &block_out);
template lila::Matrix<double>
MatrixReal<uint64_t>(BondList const &bonds, Couplings const &couplings,
                     Spinhalf<uint64_t> const &block_in,
                     Spinhalf<uint64_t> const &block_out);

template lila::Matrix<complex>
MatrixCplx<uint16_t>(BondList const &bonds, Couplings const &couplings,
                     Spinhalf<uint16_t> const &block_in,
                     Spinhalf<uint16_t> const &block_out);
template lila::Matrix<complex>
MatrixCplx<uint32_t>(BondList const &bonds, Couplings const &couplings,
                     Spinhalf<uint32_t> const &block_in,
                     Spinhalf<uint32_t> const &block_out);
template lila::Matrix<complex>
MatrixCplx<uint64_t>(BondList const &bonds, Couplings const &couplings,
                     Spinhalf<uint64_t> const &block_in,
                     Spinhalf<uint64_t> const &block_out);

} // namespace hydra
