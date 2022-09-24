#include "spinhalf_apply.h"

#include <hydra/blocks/spinhalf/terms/apply_terms_dispatch.h>
#include <hydra/blocks/spinhalf/terms/compile.h>
#include <hydra/blocks/spinhalf/terms/qns.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/operators/compiler.h>
#include <hydra/utils/logger.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
void apply(BondList const &bonds, Spinhalf<bit_t> const &block_in,
           arma::Col<coeff_t> const &vec_in, Spinhalf<bit_t> const &block_out,
           arma::Col<coeff_t> &vec_out) {

  BondList bonds_c = spinhalf::compile(bonds);
  operators::check_bonds_in_range(bonds, block_in.n_sites());

  if (block_in.n_up() != undefined_qn) {
    int n_up_out = spinhalf::nup(bonds_c) + block_in.n_up();
    if (n_up_out != block_out.n_up())
      Log.err("Incompatible n_up in Apply: {} != {}", n_up_out,
              block_out.n_up());
  }

  assert((int64_t)block_in.size() == (int64_t)vec_in.size());
  assert((int64_t)block_out.size() == (int64_t)vec_out.size());

  vec_out.zeros();

  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, coeff_t val) {
#ifdef HYDRA_ENABLE_OPENMP
    if constexpr (is_real<coeff_t>()) {
      coeff_t x = val * vec_in(idx_in);
      coeff_t *pos = vec_out.memptr();
#pragma omp atomic update
      pos[idx_out] += x;
    } else {
      complex x = val * vec_in(idx_in);
      complex *pos = vec_out.memptr();
      double *r = &reinterpret_cast<double(&)[2]>(pos[idx_out])[0];
      double *i = &reinterpret_cast<double(&)[2]>(pos[idx_out])[1];
#pragma omp atomic update
      *r += x.real();
#pragma omp atomic update
      *i += x.imag();
    }
#else
    vec_out(idx_out) += val * vec_in(idx_in);
#endif
  };

  auto const &indexing_in = block_in.indexing();
  auto const &indexing_out = block_out.indexing();
  spinhalf::apply_terms_dispatch<bit_t, coeff_t>(bonds_c, indexing_in,
                                                 indexing_out, fill);
}

template void apply<uint16_t, double>(BondList const &,
                                      Spinhalf<uint16_t> const &,
                                      arma::Col<double> const &,
                                      Spinhalf<uint16_t> const &,
                                      arma::Col<double> &);
template void apply<uint32_t, double>(BondList const &,
                                      Spinhalf<uint32_t> const &,
                                      arma::Col<double> const &,
                                      Spinhalf<uint32_t> const &,
                                      arma::Col<double> &);
template void apply<uint64_t, double>(BondList const &,
                                      Spinhalf<uint64_t> const &,
                                      arma::Col<double> const &,
                                      Spinhalf<uint64_t> const &,
                                      arma::Col<double> &);

template void apply<uint16_t, complex>(BondList const &,
                                       Spinhalf<uint16_t> const &,
                                       arma::Col<complex> const &,
                                       Spinhalf<uint16_t> const &,
                                       arma::Col<complex> &);
template void apply<uint32_t, complex>(BondList const &,
                                       Spinhalf<uint32_t> const &,
                                       arma::Col<complex> const &,
                                       Spinhalf<uint32_t> const &,
                                       arma::Col<complex> &);
template void apply<uint64_t, complex>(BondList const &,
                                       Spinhalf<uint64_t> const &,
                                       arma::Col<complex> const &,
                                       Spinhalf<uint64_t> const &,
                                       arma::Col<complex> &);

} // namespace hydra
