#include "dispatch_matrix.hpp"

#include <xdiag/algebra/fill.hpp>
#include <xdiag/basis/spinhalf/apply/dispatch.hpp>

namespace xdiag::basis::spinhalf {

template <typename coeff_t>
void dispatch_matrix(OpSum const &ops, Spinhalf const &block_in,
                     Spinhalf const &block_out, coeff_t *mat, int64_t m) try {
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_matrix(mat, idx_in, idx_out, m, val);
  };
  dispatch<coeff_t>(ops, block_in, block_out, fill);

} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
template void dispatch_matrix(OpSum const &ops, Spinhalf const &block_in,
                              Spinhalf const &block_out, double *mat,
                              int64_t m);
template void dispatch_matrix(OpSum const &ops, Spinhalf const &block_in,
                              Spinhalf const &block_out, complex *mat,
                              int64_t m);
} // namespace xdiag::basis::spinhalf
