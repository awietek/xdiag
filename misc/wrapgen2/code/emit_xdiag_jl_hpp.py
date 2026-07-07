#!/usr/bin/env python
from itertools import product
from common import *

def emit_xdiag_jl_hpp():
    total_string = """// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <utility>
#include <cstdint>

#include "jlcxx/array.hpp"
#include "jlcxx/const_array.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"
#include "jlcxx/tuple.hpp"

#include <xdiag/all.hpp>

using namespace xdiag;

// Register Mirrored types
namespace jlcxx {
"""
    # we first need to deactivate the mirrored types for plain structs

    # First, plain structs:
    with open("api.txt", "r") as fl:
        for line in [l for l in fl.readlines() if"STRUCT_DECL" in l]:
            class_name = get_class_name(line)
            total_string += "template<> struct IsMirroredType<{}> : std::false_type {{}};\n".format(class_name)

    # Second, templated structs (should be sparse matrices)
    with open("api.txt", "r") as fl:
        for l in [l for l in fl.readlines() if "CLASS_TEMPLATE" in l]:
            basename, args = class_template_basename_args(l)
            instantiations = class_template_instantiations[basename]
            combos = list(product(*instantiations.values()))
        
            for combo in combos:
                class_name, cxx_class_name = class_template_name(basename, combo)
                total_string += "template<> struct IsMirroredType<{}> : std::false_type {{}};\n".format(class_name)
    
    total_string += """}
#define JULIA_XDIAG_CALL_VOID(CMD)                                             \
  try {                                                                        \
    CMD;                                                                       \
  } catch (xdiag::Error const &e) {                                            \
    xdiag::error_trace(e);                                                     \
    throw(std::runtime_error("Error occurred in XDiag C++ core library"));     \
  }

#define JULIA_XDIAG_CALL_RETURN(CMD)                                           \
  try {                                                                        \
    return CMD;                                                                \
  } catch (xdiag::Error const &e) {                                            \
    xdiag::error_trace(e);                                                     \
    throw(std::runtime_error("Error occurred in XDiag C++ core library"));     \
  }

#define JULIA_XDIAG_CALL_RETURN_MOVE(CMD)                                      \
  try {                                                                        \
    return std::move(CMD);                                                     \
  } catch (xdiag::Error const &e) {                                            \
    xdiag::error_trace(e);                                                     \
    throw(std::runtime_error("Error occurred in XDiag C++ core library"));     \
  }

#define JULIA_XDIAG_CALL_ASSIGN(LVALUE, CMD)                                   \
  try {                                                                        \
    LVALUE = CMD;                                                              \
  } catch (xdiag::Error const &e) {                                            \
    xdiag::error_trace(e);                                                     \
    throw(std::runtime_error("Error occurred in XDiag C++ core library"));     \
  }

"""
    return total_string

print(emit_xdiag_jl_hpp())
