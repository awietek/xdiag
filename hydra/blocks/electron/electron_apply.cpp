#include "electron_apply.h"

#include <hydra/blocks/electron/terms/apply_terms_dispatch.h>
#include <hydra/blocks/electron/terms/compile.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
void apply(BondList const &bonds, Electron<bit_t> const &block_in,
           arma::Col<coeff_t> const &vec_in, Electron<bit_t> const &block_out,
           arma::Col<coeff_t> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == (idx_t)vec_in.size());
  assert(block_out.size() == (idx_t)vec_out.size());

  BondList bonds_c = electron::compile(bonds, 1e-12);

  if ((is_real<coeff_t>()) && (bonds_c.is_complex())) {
    Log.err("Error in matrix_gen: trying to create a real matrix from an "
            "intrisically complex BondList");
  }

  vec_out.zeros();
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, coeff_t val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };
  auto const &indexing_in = block_in.indexing();
  auto const &indexing_out = block_out.indexing();
  electron::apply_terms_dispatch<bit_t, coeff_t>(bonds_c, indexing_in,
                                                 indexing_out, fill);
}

template void apply<uint16_t, double>(BondList const &,
                                      Electron<uint16_t> const &,
                                      arma::Col<double> const &,
                                      Electron<uint16_t> const &,
                                      arma::Col<double> &);
template void apply<uint32_t, double>(BondList const &,
                                      Electron<uint32_t> const &,
                                      arma::Col<double> const &,
                                      Electron<uint32_t> const &,
                                      arma::Col<double> &);
template void apply<uint64_t, double>(BondList const &,
                                      Electron<uint64_t> const &,
                                      arma::Col<double> const &,
                                      Electron<uint64_t> const &,
                                      arma::Col<double> &);

template void apply<uint16_t, complex>(BondList const &,
                                       Electron<uint16_t> const &,
                                       arma::Col<complex> const &,
                                       Electron<uint16_t> const &,
                                       arma::Col<complex> &);
template void apply<uint32_t, complex>(BondList const &,
                                       Electron<uint32_t> const &,
                                       arma::Col<complex> const &,
                                       Electron<uint32_t> const &,
                                       arma::Col<complex> &);
template void apply<uint64_t, complex>(BondList const &,
                                       Electron<uint64_t> const &,
                                       arma::Col<complex> const &,
                                       Electron<uint64_t> const &,
                                       arma::Col<complex> &);

} // namespace hydra
