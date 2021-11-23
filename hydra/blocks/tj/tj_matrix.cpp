#include "tj_matrix.h"

#include <hydra/blocks/tj/terms/tj_exchange.h>
#include <hydra/blocks/tj/terms/tj_hopping.h>
#include <hydra/blocks/tj/terms/tj_ising.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
lila::Matrix<coeff_t>
MatrixGen(BondList const &bonds, Couplings const &couplings,
          tJ<bit_t> const &block_in, tJ<bit_t> const &block_out) {
  using namespace hydra::terms::tj;

  assert(block_in == block_out); // only temporary

  utils::check_operator_works_with<coeff_t>(bonds, couplings, "tj_matrix");
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  auto mat = lila::Zeros<coeff_t>(dim_out, dim_in);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, coeff_t val) {
    mat(idx_out, idx_in) += val;
  };

  auto const &indexing_in = block_in.indexing();
  // auto const &indexing_out = block_out.indexing();

  do_hopping<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  do_ising<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  do_exchange<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  return mat;
}

template lila::Matrix<double>
MatrixGen<uint16_t, double>(BondList const &bonds, Couplings const &couplings,
                            tJ<uint16_t> const &block_in,
                            tJ<uint16_t> const &block_out);
template lila::Matrix<double>
MatrixGen<uint32_t, double>(BondList const &bonds, Couplings const &couplings,
                            tJ<uint32_t> const &block_in,
                            tJ<uint32_t> const &block_out);
template lila::Matrix<double>
MatrixGen<uint64_t, double>(BondList const &bonds, Couplings const &couplings,
                            tJ<uint64_t> const &block_in,
                            tJ<uint64_t> const &block_out);

template lila::Matrix<complex>
MatrixGen<uint16_t, complex>(BondList const &bonds, Couplings const &couplings,
                             tJ<uint16_t> const &block_in,
                             tJ<uint16_t> const &block_out);
template lila::Matrix<complex>
MatrixGen<uint32_t, complex>(BondList const &bonds, Couplings const &couplings,
                             tJ<uint32_t> const &block_in,
                             tJ<uint32_t> const &block_out);
template lila::Matrix<complex>
MatrixGen<uint64_t, complex>(BondList const &bonds, Couplings const &couplings,
                             tJ<uint64_t> const &block_in,
                             tJ<uint64_t> const &block_out);

} // namespace hydra
