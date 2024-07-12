#pragma once
#ifdef XDIAG_USE_MPI
#include <xdiag/basis/spinhalf_distributed/basis_sz.h>

namespace xdiag::basis::spinhalf_distributed {

template <class bit_t, typename coeff_t>
void transpose(BasisSz<bit_t> const &basis, const coeff_t const *vec_in,
               bool reverse = true);

} // namespace xdiag::basis::spinhalf_distributed
#endif
