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
//
// Note: hermitian conjugation (hc) is intentionally NOT described here. Unlike
// the structural fields below (which are needed at Op construction/validation,
// before any algebra is known), an operator's hc behavior can carry a phase or
// even map to a different type, and is a user-facing operation. Its definition
// lives entirely in hc() (see operators/hc.cpp).
struct OpTypeInfo {
  int64_t nsites;         // required site count; undefined = no constraint
  bool site_required;     // must have a site list
  bool allow_overlapping; // repeated sites are allowed
  bool matrix_required;   // must carry an arma matrix
  bool is_real;           // matrix is always real
};

// clang-format off
inline const std::map<std::string, OpTypeInfo> op_registry = {
//  type             nsites      site_req  overlap  mat_req  is_real
  {"Id",            {undefined,  false,    false,   false,   true   }},
  {"SdotS",         {2,          true,     true,    false,   true   }},
  {"Exchange",      {2,          true,     true,    false,   true   }},
  {"ExchangeAsym",  {2,          true,     true,    false,   true   }},
  {"SzSz",          {2,          true,     true,    false,   true   }},
  {"Sz",            {1,          true,     false,   false,   true   }},
  {"S+",            {1,          true,     false,   false,   true   }},
  {"S-",            {1,          true,     false,   false,   true   }},
  {"ScalarChirality",{3,         true,     false,   false,   false  }},
  {"Hop",           {2,          true,     false,   false,   true   }},
  {"Hopup",         {2,          true,     false,   false,   true   }},
  {"Hopdn",         {2,          true,     false,   false,   true   }},
  {"Cdagup",        {1,          true,     false,   false,   true   }},
  {"Cup",           {1,          true,     false,   false,   true   }},
  {"Cdagdn",        {1,          true,     false,   false,   true   }},
  {"Cdn",           {1,          true,     false,   false,   true   }},
  {"HubbardU",      {undefined,  false,    false,   false,   true   }},
  {"Ntot",          {1,          true,     false,   false,   true   }},
  {"Nup",           {1,          true,     false,   false,   true   }},
  {"Ndn",           {1,          true,     false,   false,   true   }},
  {"Nupdn",         {1,          true,     false,   false,   true   }},
  {"NtotNtot",      {2,          true,     false,   false,   true   }},
  {"NupdnNupdn",    {2,          true,     false,   false,   true   }},
  {"tJSzSz",        {2,          true,     false,   false,   true   }},
  {"tJSdotS",       {2,          true,     false,   false,   true   }},
  {"Matrix",        {undefined,  true,     true,    true,    false  }},
  {"NupNdn",        {2,          true,     false,   false,   true   }},
  {"NupNup",        {2,          true,     false,   false,   true   }},
  {"NdnNdn",        {2,          true,     false,   false,   true   }},
  {"NdnNup",        {2,          true,     false,   false,   true   }},
};
// clang-format on

// Query functions
bool is_known_type(std::string const &type);
bool is_real_type(std::string const &type);
int64_t nsites_of_type(std::string const &type);
OpTypeInfo const &info_of_type(std::string const &type);
std::string known_types_string();

} // namespace xdiag
