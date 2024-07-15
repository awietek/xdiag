#include "tj_distributed_apply.hpp"

#include <xdiag/algebra/fill.hpp>
#include <xdiag/blocks/tj/compile.hpp>
#include <xdiag/blocks/tj_distributed/dispatch.hpp>

namespace xdiag {

template <typename coeff_t>
void apply(BondList const &bonds, tJDistributed const &block_in,
           arma::Col<coeff_t> const &vec_in, tJDistributed const &block_out,
           arma::Col<coeff_t> &vec_out, double zero_precision) try {
  vec_out.zeros();
  BondList bondsc = tj::compile(bonds, block_in.n_sites(), zero_precision);
  tj_distributed::dispatch<coeff_t>(bondsc, block_in, vec_in, block_out,
                                    vec_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply<double>(BondList const &, tJDistributed const &,
                            arma::Col<double> const &, tJDistributed const &,
                            arma::Col<double> &, double);

template void apply<complex>(BondList const &, tJDistributed const &,
                             arma::Col<complex> const &, tJDistributed const &,
                             arma::Col<complex> &, double);

} // namespace xdiag
