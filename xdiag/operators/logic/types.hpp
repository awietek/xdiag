// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <limits>
#include <map>
#include <string>

namespace xdiag {

constexpr int64_t undefined = std::numeric_limits<int64_t>::min();

// Per-type static properties
struct OpTypeInfo {
  int64_t     nsites;            // required site count; undefined = no constraint
  bool        site_required;     // must have a site list
  bool        allow_overlapping; // repeated sites are allowed
  bool        matrix_required;   // must carry an arma matrix
  bool        is_real;           // matrix is always real
  std::string hc_partner;        // same name = self-adjoint; different name = partner
                                 // for "Matrix": also conjugate-transposes the matrix
};

// clang-format off
inline const std::map<std::string, OpTypeInfo> op_registry = {
//  type             nsites      site_req  overlap  mat_req  is_real  hc_partner
  {"Id",            {undefined,  false,    false,   false,   true,    "Id"           }},
  {"SdotS",         {2,          true,     true,    false,   true,    "SdotS"        }},
  {"Exchange",      {2,          true,     true,    false,   true,    "Exchange"     }},
  {"SzSz",          {2,          true,     true,    false,   true,    "SzSz"         }},
  {"Sz",            {1,          true,     false,   false,   true,    "Sz"           }},
  {"S+",            {1,          true,     false,   false,   true,    "S-"           }},
  {"S-",            {1,          true,     false,   false,   true,    "S+"           }},
  {"ScalarChirality",{3,         true,     false,   false,   false,   "ScalarChirality"}},
  {"Hop",           {2,          true,     false,   false,   true,    "Hop"          }},
  {"Hopup",         {2,          true,     false,   false,   true,    "Hopup"        }},
  {"Hopdn",         {2,          true,     false,   false,   true,    "Hopdn"        }},
  {"Cdagup",        {1,          true,     false,   false,   true,    "Cup"          }},
  {"Cup",           {1,          true,     false,   false,   true,    "Cdagup"       }},
  {"Cdagdn",        {1,          true,     false,   false,   true,    "Cdn"          }},
  {"Cdn",           {1,          true,     false,   false,   true,    "Cdagdn"       }},
  {"HubbardU",      {undefined,  false,    false,   false,   true,    "HubbardU"     }},
  {"Ntot",          {1,          true,     false,   false,   true,    "Ntot"         }},
  {"Nup",           {1,          true,     false,   false,   true,    "Nup"          }},
  {"Ndn",           {1,          true,     false,   false,   true,    "Ndn"          }},
  {"Nupdn",         {1,          true,     false,   false,   true,    "Nupdn"        }},
  {"NtotNtot",      {2,          true,     false,   false,   true,    "NtotNtot"     }},
  {"NupdnNupdn",    {2,          true,     false,   false,   true,    "NupdnNupdn"   }},
  {"tJSzSz",        {2,          true,     false,   false,   true,    "tJSzSz"       }},
  {"tJSdotS",       {2,          true,     false,   false,   true,    "tJSdotS"      }},
  {"Matrix",        {undefined,  true,     true,    true,    false,   "Matrix"       }},
  {"NupNdn",        {2,          true,     false,   false,   true,    "NupNdn"       }},
  {"NupNup",        {2,          true,     false,   false,   true,    "NupNup"       }},
  {"NdnNdn",        {2,          true,     false,   false,   true,    "NdnNdn"       }},
  {"NdnNup",        {2,          true,     false,   false,   true,    "NdnNup"       }},
};
// clang-format on

// Query functions
bool is_known_type(std::string const &type);
bool is_real_type(std::string const &type);
int64_t nsites_of_type(std::string const &type);
OpTypeInfo const &info_of_type(std::string const &type);
std::string known_types_string();

} // namespace xdiag
