#include "tj_apply.h"

#include <hydra/blocks/tj/terms/apply_terms_dispatch.h>
#include <hydra/blocks/tj/terms/compile.h>
#include <hydra/operators/compiler.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
void apply(BondList const &bonds, tJ<bit_t> const &block_in,
           arma::Col<coeff_t> const &vec_in, tJ<bit_t> const &block_out,
           arma::Col<coeff_t> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert((idx_t)block_in.size() == (idx_t)vec_in.size());
  assert((idx_t)block_out.size() == (idx_t)vec_out.size());

  BondList bonds_c = tj::compile(bonds, 1e-12);
  operators::check_bonds_in_range(bonds, block_in.n_sites());

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
  tj::apply_terms_dispatch<bit_t, coeff_t>(bonds_c, indexing_in, indexing_out,
                                           fill);
}

template void apply<uint16_t, double>(BondList const &, tJ<uint16_t> const &,
                                      arma::Col<double> const &,
                                      tJ<uint16_t> const &,
                                      arma::Col<double> &);
template void apply<uint32_t, double>(BondList const &, tJ<uint32_t> const &,
                                      arma::Col<double> const &,
                                      tJ<uint32_t> const &,
                                      arma::Col<double> &);
template void apply<uint64_t, double>(BondList const &, tJ<uint64_t> const &,
                                      arma::Col<double> const &,
                                      tJ<uint64_t> const &,
                                      arma::Col<double> &);

template void apply<uint16_t, complex>(BondList const &, tJ<uint16_t> const &,
                                       arma::Col<complex> const &,
                                       tJ<uint16_t> const &,
                                       arma::Col<complex> &);
template void apply<uint32_t, complex>(BondList const &, tJ<uint32_t> const &,
                                       arma::Col<complex> const &,
                                       tJ<uint32_t> const &,
                                       arma::Col<complex> &);
template void apply<uint64_t, complex>(BondList const &, tJ<uint64_t> const &,
                                       arma::Col<complex> const &,
                                       tJ<uint64_t> const &,
                                       arma::Col<complex> &);

} // namespace hydra
