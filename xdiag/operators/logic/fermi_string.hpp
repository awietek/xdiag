// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <tuple>

#include <xdiag/operators/op.hpp>

namespace xdiag::operators {

bool is_fermi_string(Op const &op);
int64_t check_fermi_string(std::string type);

enum FermiOperator { cdagup, cdagdn, cup, cdn };

template <typename bit_t> class FermiString {
public:
  FermiString() = default;
  FermiString(Op const &op);
  bool non_zero_term(bit_t up, bit_t dn) const;
  std::tuple<bit_t, bit_t, double> up_dn_coeff(bit_t up, bit_t dn) const;

private:
  int64_t size_;
  std::vector<FermiOp> operator_;
  std::vector<bit_t> set_mask_;
  std::vector<bit_t> fermi_mask_;
};

} // namespace xdiag::operators
