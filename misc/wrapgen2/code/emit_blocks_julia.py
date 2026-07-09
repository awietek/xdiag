#!/usr/bin/env python
import sys
from common import *

def emit_blocks_julia():
    total_string = """# SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0
# automatically generated using emit_blocks_julia.py

abstract type Block end
    
"""
    module_lines = get_module_lines("api.txt", "blocks")
    for block in blocks:
        total_string += """#####################################
# {} block type
""".format(block)
        
        total_string += """struct {0} <: Block
    cxx_object::typ_{0}
end
convert(::Type{{T}}, block::typ_{0}) where T <: {0} = {0}(block)
        
""".format(block)

        total_string += "# Constructors\n"
        for line in module_lines:
            if "CONSTRUCTOR" in line:
                class_name = get_parent_class_name(line)
                if class_name != block:
                    continue

                args_list, args_raw_list, defaults_list = get_args_lists(line)
                args, args_raw = get_args(line)                
                # the block constructors with RepresentationSet
                # will not be exposed to Julia
                if "RepresentationSet" in args:
                    continue
                total_string += "{0}({1}) = {0}(con_{0}({2}))\n".format(
                    block, args_to_julia(args_list, defaults_list), args_raw)

        total_string += "\n"

        total_string += "# Methods\n"
        for line in module_lines:
            if "CXX_METHOD" in line:
                class_name = get_parent_class_name(line)
                if class_name != block:
                    continue

                method_name, return_type = get_name_return_type(line, "CXX_METHOD")
                if method_name in ignored_methods or return_type in ignored_types:
                    continue
                
                args_list, args_raw_list, defaults_list = get_args_lists(line)
                args, args_raw = get_args(line)
                full_method_name = method_name_with_base(method_name)
    
                if args == "":
                    total_string += "{}(self::{})::{} = mth_{}(self.cxx_object)\n".format(
                        full_method_name, class_name,
                        return_type_to_julia(return_type), method_name)
                else:
                    total_string += "{}(self::{}, {})::{} = mth_{}(self.cxx_object, {})\n".format(
                        full_method_name, class_name,
                        args_to_julia(args_list, defaults_list),
                        return_type_to_julia(return_type), method_name,
                        args_raw)
        total_string += "\n"

        
#         total_string += """# Iterators
# function Base.iterate(block::{0})
#     if size(block) > 0
#         b = _begin(block.cxx_object)
#         return ProductState(_deref(b)), b
#     else
#         return nothing
#     end
# end

# function Base.iterate(block::{0}, state)
#     _incr(state)
#     if state != _end(block.cxx_object)
#         return ProductState(_deref(state)), state
#     else
#         return nothing
#     end
# end
# """
        
        total_string += """# Output
to_string(block::{0})::String = fun_to_string(block.cxx_object)
Base.show(io::IO, block::{0}) = print(io, fun_to_string(block.cxx_object))
        
""".format(block)
        
    return total_string

# emit_blocks_julia()
print(emit_blocks_julia())

