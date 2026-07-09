#!/usr/bin/env python
import sys
from common import *
        
def emit_types():
    ts = """# SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0
# automatically generated using emit_types_julia.py
"""              

    with open("api.txt", 'r') as fl:
        for line in fl.readlines():

            if is_deleted(line):
                continue
        
            if "CLASS_DECL" in line or "STRUCT_DECL" in line:
                class_name = get_class_name(line)
                ts += """
struct {0}
    cxx_object::typ_{0}
end
Base.convert(::Type{{{0}}}, obj::typ_{0}Allocated) = {0}(obj)   
""".format(class_name)

    return ts
                
print(emit_types())
