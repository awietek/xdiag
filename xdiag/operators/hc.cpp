// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "hc.hpp"

#include <map>
#include <string>
#include <utility>

#include <xdiag/operators/valid.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag {

// Hermitian-conjugate partner of each Op type: hc(type{sites}) is
//
//     phase * partner{sites}
//
// where {phase, partner} is looked up below. Self-adjoint types map to
// themselves with phase +1; raising/lowering and creation/annihilation types
// map to their counterpart; ExchangeAsym is anti-hermitian and carries -1.
//
// Types absent from this table have no well-defined hermitian conjugate, and
// hc() throws on them. The "Matrix" type is handled separately (its stored
// matrix is conjugate-transposed), so it is intentionally not listed here.
//
// The phase is real for every operator type here (it is only ever +1 or -1); a
// real phase also keeps hc of a real operator real. Should a type ever require
// a genuinely complex phase, widen this to complex.
//
// clang-format off
static const std::map<std::string, std::pair<double, std::string>> hc_partners = {
//  type                phase   partner
  {"Id",              {1.0,  "Id"             }},
  {"SdotS",           {1.0,  "SdotS"          }},
  {"Exchange",        {1.0,  "Exchange"       }},
  {"ExchangeAsym",    {-1.0, "ExchangeAsym"   }}, // anti-hermitian
  {"SzSz",            {1.0,  "SzSz"           }},
  {"Sz",              {1.0,  "Sz"             }},
  {"S+",              {1.0,  "S-"             }},
  {"S-",              {1.0,  "S+"             }},
  {"ScalarChirality", {1.0,  "ScalarChirality"}},
  {"Hop",             {1.0,  "Hop"            }},
  {"Hopup",           {1.0,  "Hopup"          }},
  {"Hopdn",           {1.0,  "Hopdn"          }},
  {"Cdagup",          {1.0,  "Cup"            }},
  {"Cup",             {1.0,  "Cdagup"         }},
  {"Cdagdn",          {1.0,  "Cdn"            }},
  {"Cdn",             {1.0,  "Cdagdn"         }},
  {"HubbardU",        {1.0,  "HubbardU"       }},
  {"Ntot",            {1.0,  "Ntot"           }},
  {"Nup",             {1.0,  "Nup"            }},
  {"Ndn",             {1.0,  "Ndn"            }},
  {"Nupdn",           {1.0,  "Nupdn"          }},
  {"NtotNtot",        {1.0,  "NtotNtot"       }},
  {"NupdnNupdn",      {1.0,  "NupdnNupdn"     }},
  {"tJSzSz",          {1.0,  "tJSzSz"         }},
  {"tJSdotS",         {1.0,  "tJSdotS"        }},
  {"NupNdn",          {1.0,  "NupNdn"         }},
  {"NupNup",          {1.0,  "NupNup"         }},
  {"NdnNdn",          {1.0,  "NdnNdn"         }},
  {"NdnNup",          {1.0,  "NdnNup"         }},
};
// clang-format on

OpSum hc(Op const &op) try {
  operators::check_valid(op);

  // "Matrix" type: same type and sites, but the matrix is conjugate-transposed.
  if (op.type() == "Matrix") {
    return OpSum(Op("Matrix", op.sites(), op.matrix().hc()));
  }

  auto it = hc_partners.find(op.type());
  if (it == hc_partners.end()) {
    XDIAG_THROW(fmt::format(
        "Cannot form the hermitian conjugate (hc) of Op type \"{}\": no "
        "hermitian-conjugate partner is defined for it.",
        op.type()));
  }

  auto const &[phase, partner] = it->second;
  if (op.hassites()) {
    return OpSum(phase, Op(partner, op.sites()));
  } else {
    return OpSum(phase, Op(partner));
  }
}
XDIAG_CATCH

OpSum hc(Monomial const &mono) try {
  // hc(A_1 * A_2 * ... * A_n) = hc(A_n) * ... * hc(A_2) * hc(A_1)
  OpSum mono_hc(Monomial{}); // multiplicative identity (empty product = 1)
  for (int64_t i = mono.size() - 1; i >= 0; --i) {
    mono_hc = mono_hc * hc(mono[i]);
  }
  return mono_hc;
}
XDIAG_CATCH

OpSum hc(OpSum const &ops) try {
  OpSum ops_hc;
  for (auto const &[coeff, mono] : ops.plain()) {
    // hc(coeff * mono) = conj(coeff) * hc(mono)
    ops_hc += conj(coeff.scalar()) * hc(mono);
  }
  return ops_hc;
}
XDIAG_CATCH

} // namespace xdiag
