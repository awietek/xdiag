#!/usr/bin/env python
import os
from common import *

def emit_xdiag_jl():
    ts = """# SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0
# automatically generated using create_wrapper_julia.py

module XDiag
using CxxWrap
using LinearAlgebra
using XDiag_jll
    
"""
    ts += "import Base: {}\n".format(", ".join(julia_base_imports))
    ts += "import LinearAlgebra: {}\n\n".format(", ".join(julia_linalg_imports))
    ts += "export {}\n".format(", ".join(julia_exports))
    
    
    ts += """
@wrapmodule(XDiag_jll.get_libxdiagjl_path)

include("armadillo.jl")
include("types.jl")
"""
    modules = [os.path.basename(f.path) for f in os.scandir(os.path.join(project_root, "xdiag")) if f.is_dir()]
    for mod in modules:
        ts += "include(\"{}.jl\")\n".format(mod)

    ts +="""
function __init__()
    @initcxx
end

end
"""
    return ts
    
print(emit_xdiag_jl())


