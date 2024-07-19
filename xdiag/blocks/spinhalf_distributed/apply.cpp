#include "apply.hpp"

#include <xdiag/basis/spinhalf_distributed/apply.hpp>
#include <xdiag/blocks/spinhalf/compile.hpp>

namespace xdiag {
template <typename coeff_t>
void apply(OpSum const &ops, SpinhalfDistributed const &block_in,
           arma::Col<coeff_t> const &vec_in,
           SpinhalfDistributed const &block_out, arma::Col<coeff_t> &vec_out,
           double zero_precision) try {
  Log("Arrived block");
  vec_out.zeros();
  OpSum opsc = spinhalf::compile(ops, block_in.n_sites(), zero_precision);

  std::visit(
      [&](auto &&basis_in, auto &&basis_out) {
        using basis_in_t = typename std::decay<decltype(basis_in)>::type;
        using basis_out_t = typename std::decay<decltype(basis_out)>::type;
        if constexpr (std::is_same<basis_in_t, basis_out_t>::value) {
          basis::spinhalf_distributed::apply(ops, basis_in, vec_in, basis_out,
                                             vec_out);
        } else {
          XDIAG_THROW(
              "Invalid combination of bases for \"SpinhalfDistributed\" block.")
        }
      },
      block_in.basis(), block_out.basis());

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply<double>(OpSum const &, SpinhalfDistributed const &,
                            arma::Col<double> const &,
                            SpinhalfDistributed const &, arma::Col<double> &,
                            double);

template void apply<complex>(OpSum const &, SpinhalfDistributed const &,
                             arma::Col<complex> const &,
                             SpinhalfDistributed const &, arma::Col<complex> &,
                             double);

} // namespace xdiag
