#!/usr/bin/env python
import sys
from common import *
        
def emit_submodule(module):
    ts = """# SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0
# automatically generated using emit_submodule_julia.py
"""              
    module_lines = get_module_lines("api.txt", module)
    method_names = []
    for line in module_lines:
        if "CLASS_DECL" in line:
            class_name = get_class_name(line)
            ts += "\n# {}\n".format(class_name)
        
        if "CONSTRUCTOR" in line:
            class_name = get_parent_class_name(line)

            if class_name in ignored_classes:
                continue

            if contains_invalid_type(line):
                continue

            args_to_julia, args_to_cpp = args_to_julia_to_cpp(line)
            ts += "{0}({1}) = {0}(con_{0}({2}))\n".format(
                class_name, args_to_julia, args_to_cpp)

        elif "CXX_METHOD" in line:
            class_name = get_parent_class_name(line)
            method_name, return_type = get_name_return_type(line, "CXX_METHOD")
            if return_type not in julia_types.keys():
                continue

            if method_name in ignored_methods or return_type in ignored_types:
                continue

            if contains_invalid_type(line):
                continue

            full_method_name = method_name_with_base(method_name)
            args_to_julia, args_to_cpp = args_to_julia_to_cpp(line)
            if args_to_julia == "":
                ts += "{}(obj::{})::{} = mth_{}(obj.cxx_object)\n".format(
                    full_method_name, class_name, julia_types[return_type], method_name)
            else:
                ts += "{}(obj::{}, {})::{} = mth_{}(obj.cxx_object, {})\n".format(
                    full_method_name, class_name, args_to_julia,
                    method_name, julia_types[return_type], args_to_cpp)

            method_names.append(method_name)

    # we first want all method names to be classes and method_names
    # to be collected before we move to free functions
    ts += "\n# Free functions \n"
    for line in module_lines:

        if "FUNCTION_DECL" in line:
            function_name, return_type = get_name_return_type(line, "FUNCTION_DECL")
            if return_type not in julia_types.keys():
                continue 

            if function_name in ignored_methods or return_type in ignored_types:
                continue

            if contains_invalid_type(line):
                continue

            # Avoid collisions between free functions and methods in C++
            if function_name in method_names:
                continue

            full_function_name = method_name_with_base(function_name)
            args_to_julia, args_to_cpp = args_to_julia_to_cpp(line)
            ts += "{}({})::{} = fun_{}({})\n".format(
                    full_function_name, args_to_julia,
                    julia_types[return_type], function_name, args_to_cpp)
            if function_name == "to_string":
                ts += "Base.show(io::IO, {0}) = print(io, fun_to_string({1}))\n".format(args_to_julia, args_to_cpp)


    # Finally, operator overloads
    ts += "\n# operator overloads \n"
    for line in module_lines:
        if "CXX_METHOD" in line:
            class_name = get_parent_class_name(line)
            method_name, return_type = get_name_return_type(line, "CXX_METHOD")

            # operator[] get's special treatment
            if(method_name == "operator[]"):
                # Ignore all weird return types except "Scalar"
                if return_type in ignored_types and not "Scalar" in return_type:
                    continue                    
                
                args, args_raw = get_args(line)
                args_to_julia, args_to_cpp = args_to_julia_to_cpp(line)
                if (is_const_qualified(line)):
                    if "Scalar" not in return_type:    # read access
                        ts += "Base.getindex(self::{}, {}) = op_getindex(self.cxx_object, {})\n".format(class_name, args_to_julia, args_to_cpp)
                else:                             # write access
                    raw_return_type = decay_type(return_type)
                    if raw_return_type == "Scalar":
                        types = ["double", "complex"]
                    else:
                        types = [raw_return_type]
                    for typ in types:
                        ts += "Base.setindex!(self::{}, val::{}, {}) = op_setindex(self.cxx_object, {}, val)\n".format(class_name, julia_types[typ], args_to_julia, args_to_cpp)
                continue
            
            if return_type not in julia_types.keys() or\
               return_type in ignored_types:
                continue

            if method_name not in operator_aux_julia_name.keys():
                continue

            if contains_invalid_type(line):
                continue

            aux_method_name, julia_method_name = operator_aux_julia_name[method_name]
            args_to_julia, args_to_cpp = args_to_julia_to_cpp(line)
            ts += "{}(self::{}, {})::{} = op_{}(self.cxx_object, {})\n".format(
                julia_method_name, class_name, args_to_julia,
                julia_types[return_type], aux_method_name, args_to_cpp)

        elif "FUNCTION_DECL" in line:
            method_name, return_type = get_name_return_type(line, "FUNCTION_DECL")
            if return_type not in julia_types.keys() or\
               return_type in ignored_types:
                continue

            if method_name not in operator_aux_julia_name.keys():
                continue

            if contains_invalid_type(line):
                continue
            
            aux_method_name, julia_method_name = operator_aux_julia_name[method_name]
            args_to_julia, args_to_cpp = args_to_julia_to_cpp(line)

            ts += "{}({})::{} = op_{}({})\n".format(
                julia_method_name, args_to_julia, julia_types[return_type],
                aux_method_name, args_to_cpp)
            
    return ts

module = str(sys.argv[1])
# emit_submodule(module)
print(emit_submodule(module))
