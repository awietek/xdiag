// SPDX-FileCopyrightText: 2025 
//
// SPDX-License-Identifier: Apache-2.0 

#pragma once 
#include <xdiag/common.hpp>
#include <xdiag/operators/coupling.hpp> 

namespace xdiag::basis {

template <typename coeff_t, class basis_t>
void apply_identity_distributed(Coupling const &cpl, basis_t const &basis, 
                                const coeff_t *vec_in, coeff_t *vec_out) {
    using bit_t = typename basis_t::bit_t;
    coeff_t s = cpl.scalar().as<coeff_t>();

    int64_t idx = 0;
    for (bit_t up : basis.my_ups()) {
        for (bit_t dn : basis.all_dns()) {
            vec_out[idx] += s * vec_in[idx];
            ++idx;
        }
    }
}
} // namespace xdiag::basis 
