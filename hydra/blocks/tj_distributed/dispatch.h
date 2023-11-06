#pragma once
#ifdef HYDRA_USE_MPI

#include <hydra/common.h>

#include <hydra/basis/basis.h>
#include <hydra/blocks/tj_distributed/terms/apply_terms.h>

namespace hydra::tj_distributed {

template <typename coeff_t, class fill_f>
inline void dispatch(BondList const &bonds, tJDistributed const &block_in,
                     tJDistributed const &block_out, fill_f &&fill) try {
  using namespace basis::tj_distributed;

  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(
      overload{// uint16_t
               [&](BasisNp<uint16_t> const &idx_in,
                   BasisNp<uint16_t> const &idx_out) {
                 apply_terms<uint16_t, coeff_t>(bonds, idx_in, idx_out, fill);
               },
               [&](BasisNp<uint32_t> const &idx_in,
                   BasisNp<uint32_t> const &idx_out) {
                 apply_terms<uint32_t, coeff_t>(bonds, idx_in, idx_out, fill);
               },
               [&](BasisNp<uint64_t> const &idx_in,
                   BasisNp<uint64_t> const &idx_out) {
                 apply_terms<uint64_t, coeff_t>(bonds, idx_in, idx_out, fill);
               },
               [&](auto const &idx_in, auto const &idx_out) {
                 HydraThrow(std::logic_error,
                            "Invalid basis or combination of bases");
                 (void)idx_in;
                 (void)idx_out;
               }},
      basis_in, basis_out);
} catch (...) {
  HydraRethrow("Unable to apply terms on tJ block");
}

} // namespace hydra::tj_distributed

#endif
