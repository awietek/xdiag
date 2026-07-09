project_root = "/Users/awietek/Research/Software/xdiag"

class_template_instantiations = \
    {"COOMatrix": {"idx_t": ["int32_t", "int64_t"],
                   "coeff_t": ["double", "complex"]},
     "CSRMatrix": {"idx_t": ["int32_t", "int64_t"],
                   "coeff_t": ["double", "complex"]},
     "CSCMatrix": {"idx_t": ["int32_t", "int64_t"],
                   "coeff_t": ["double", "complex"]}}

julia_types = {"bool" : "Bool",
               "int64_t": "Int64",
               "double": "Float64",
               "complex": "ComplexF64",
               "std::string": "String",
               "arma::Col<int64_t>": "Vector{Int64}",
               "arma::Mat<int64_t>": "Matrix{Int64}",
               "arma::vec": "Vector{Float64}",
               "arma::cx_vec": "Vector{ComplexF64}",
               "arma::mat": "Matrix{Float64}",
               "arma::cx_mat": "Matrix{ComplexF64}",
               "Permutation": "Permutation",
               "PermutationGroup": "PermutationGroup",
               "Representation": "Representation",
               "Op": "Op",
               "Monomial": "Monomial",
               "OpSum": "OpSum",
               "ProductState": "ProductState"
               }

blocks = ["Fermion", "Boson", "Spinhalf", "Electron", "tJ"]

julia_base_imports = ["getindex", "setindex!", "size", "isreal", "convert",
                      "show", "real", "imag", "push!", "iterate", "length",
                      "isapprox", "inv", "conj", "abs", "zero", "adjoint",
                      "+", "-", "*", "/"]

julia_linalg_imports = ["dot", "norm"]

julia_exports = ["say_hello", "print_version", "set_verbosity",
                 "OpSum", "Op", "Monomial", "plain"] + blocks

operator_aux_julia_name = {"operator+": ("add", "Base.:+"),
                           "operator-": ("sub", "Base.:-"),
                           "operator*": ("mul", "Base.:*"),
                           "operator/": ("div", "Base.:/"),
                           "operator==": ("isequal", "Base.:(==)"),
                           "operator++": ("increment", "increment")}
                           
ignored_methods = ["fill", "what",
                   "operator+", "operator-", "operator*", "operator/",
                   "operator+=", "operator-=", "operator*=", "operator/=",
                   "operator==", "operator++", "operator<<", "operator[]",
                   "operator!=", "operator<",
                   "begin", "end",
                   "irreps", # in blocks: would return internal RepresenationSet
                   "basis", # in blocks: would return std::shared_ptr
                   "table", # in FileToml: would return the raw toml table
                   "terms", "params", # part of OpSum: for iteration
                   "block",           # part of State
                   "rightnow"   # timing not necessary for Julia
                   ]
ignored_classes = ["FileTomlHandler", "GPWF"]
ignored_types = ["std::istream", "std::initializer_list",
                 "Matrix", "Matrix const &", "Matrix const&",
                 "Vector", "Vector const &", "Vector const&",
                 "Scalar", "Scalar const &", "Scalar const&",
                 "Coeff", "Coeff const &", "Coeff const&",
                 "Term", "Term const &", "Term const&",
                 "Error", "FileTomlHandler"]

def decay_type(typ):
    return typ.replace("const", "").replace("&", "").strip()

def is_deleted(line):
    return line.split("=")[-1].strip() == "delete"

def all_xdiag_types():
    class_names = []
    with open("api.txt", "r") as fl:
        for line in [l for l in fl.readlines() if "CLASS_DECL" in l or "STRUCT_DECL" in l]:
            class_names.append(get_class_name(line))
    return class_names

def is_const_qualified(line):
    return line.split(")")[-1].strip() == "const"

def args_to_julia(args_list, defaults_list):
    julia_args = []
    for arg, default in zip(args_list, defaults_list):
        clean_arg = arg.replace("const &", "").replace("const&", "")
        spt = clean_arg.split()
        if spt[0] in julia_types:
            julia_arg = "{}::{}".format(spt[1], julia_types[spt[0]])
        else:
            julia_arg = "{}::{}".format(spt[1], spt[0])
        if default != None:
            julia_arg += " = {}".format(default)
        julia_args.append(julia_arg)
    return ", ".join(julia_args)

def args_to_julia_to_cpp(line):
    args_list, args_raw_list, types_list, defaults_list = get_args_types_defaults(line)
    args_to_julia = []
    args_to_cpp = []
    for arg, typ, default in zip(args_raw_list, types_list, defaults_list):
        typ = typ.replace("const", "").replace("&", "").strip()
        
        if typ in julia_types.keys():
            julia_typ = julia_types[typ]
        else:
            julia_typ = typ
        if default == None:
            args_to_julia.append("{}::{}".format(arg, julia_typ))
        else:
            args_to_julia.append("{}::{} = {}".format(arg, julia_typ,
                                                      default))
                    
        if "arma" in typ:
            args_to_cpp.append("to_armadillo({})".format(arg))
        elif typ in all_xdiag_types():
            args_to_cpp.append("{}.cxx_object".format(arg))
        else:
            args_to_cpp.append("{}".format(arg))

    args_to_julia = ", ".join(args_to_julia)
    args_to_cpp = ", ".join(args_to_cpp)
    return args_to_julia, args_to_cpp

def contains_invalid_type(line):
    args_list, args_raw_list, types_list, defaults_list = get_args_types_defaults(line)
    for arg, typ, default in zip(args_raw_list, types_list, defaults_list):
        typ = typ.replace("const", "").replace("&", "").strip()
        if typ not in julia_types.keys() and typ not in all_xdiag_types():
            return True
    return False
                
def return_type_to_julia(typ):
    if typ in julia_types:
        return julia_types[typ]
    else:
        return typ

def method_name_with_base(method_name):
    if method_name in julia_base_imports:
        return "Base.{}".format(method_name)
    elif method_name in julia_linalg_imports:
        return "LinearAlgebra.{}".format(method_name)
    else:
        return method_name
    
def get_module_lines(filename, module):
    module_lines = []
    is_module_line = False
    total_header = "=== xdiag"
    module_header = "=== xdiag/{}".format(module)

    with open("api.txt", "r") as fl:
        for line in fl.readlines():
            if is_module_line:
                if total_header in line and module_header not in line:
                    is_module_line = False
                    break
                module_lines.append(str(line).strip())
                
            else:   # !is_module_line
                if module_header in line:
                    is_module_line = True
    return module_lines

# If some argument list contains a "Block", then we cannot expose it with the
# Block interface, since CxxWrap can't handle the std::variant
# Instead each Block gets replaced by the concrete block types
def get_substituted_args(args):
    for t in ignored_types:
        if t in args:
            return []
        
    if "Block" in args:
        return [b.join(args.split("Block")) for b in blocks]
    else:
        return [args]

def get_parent_class_name(line):
    return str(line.split("[")[1].split("]")[0]).strip()

def get_args_lists(line):
    args_with_defaults = line.split("(")[1].split(")")[0].split(",")
    args_list = []
    args_raw_list = []
    defaults_list = []
    for a in args_with_defaults:
        if a == "":
            continue
        
        spt = a.split("=")
        arg = spt[0].strip()
        arg_raw = arg.split(" ")[-1].replace("&", "").replace("*", "")
        args_list.append(arg)
        args_raw_list.append(arg_raw)
        if len(spt) == 1:
            defaults_list.append(None)
        else:
            defaults_list.append(a.split("=")[1].strip())
    return args_list, args_raw_list, defaults_list


def get_args_types_defaults(line):
    args_with_defaults = line.split("(")[1].split(")")[0].split(",")
    args_list = []
    args_raw_list = []
    types_list = []
    defaults_list = []
    for a in args_with_defaults:
        if a == "":
            continue
        
        spt = a.split("=")
        arg = spt[0].strip()
        typ = " ".join(arg.split(" ")[:-1])
        arg_raw = arg.split(" ")[-1]
        if arg_raw[0] == "&":
            arg_raw = arg_raw[1:]
            typ = typ + "&"
        if arg_raw[0] == "*":
            arg_raw = arg_raw[1:]
            typ = typ + "*"

        args_list.append(arg)
        args_raw_list.append(arg_raw)
        types_list.append(typ)
        if len(spt) == 1:
            defaults_list.append(None)
        else:
            defaults_list.append(a.split("=")[1].strip())
    return args_list, args_raw_list, types_list, defaults_list


def get_args(line):
    args_list, args_raw_list, _ = get_args_lists(line)
    args = ", ".join(args_list)
    args_raw = ", ".join(args_raw_list)
    return args, args_raw

def get_name_return_type(line, apitype):
    method_name = line.split("(")[0].split(" ")[-1].strip()
    return_type = " ".join(line.split(apitype)[1].split("(")[0].split(" ")[:-1]).replace("constexpr", "").replace("XDIAG_API", "").strip()

    # correct for reference
    if method_name[0] == "&":
        method_name = method_name[1:]
        return_type += "&"

    # correct for pointer
    if method_name[0] == "*":
        method_name = method_name[1:]
        return_type += "*"
    return method_name, return_type
        

def class_template_basename_args(line):
    head = line.split("{")[0].strip()
    basename = head.split(" ")[-1]
    templateargs_str = head.split("<")[1].split(">")[0].strip()
    templateargs = [str(a.strip()) for a in templateargs_str.replace("typename", " ").split(",")]
    return basename, templateargs
    
def class_template_name(basename, combo):
    cxx_class_name = "_".join(["typ", basename] + list(combo))
    tplargs = ",".join(combo)
    class_name = basename + "<" + tplargs + ">"
    return class_name, cxx_class_name

def get_class_name(line):
    # class name comes directly after XDIAG_API
    spt = line.split(" ")
    idx = spt.index("XDIAG_API")
    assert str(spt[idx+2]) == "{"
    return spt[idx+1]

def get_struct_name_fields(line):
    # class name comes directly after XDIAG_API
    spt = line.split(" ")
    idx = spt.index("XDIAG_API")
    name = spt[idx+1]
    assert str(spt[idx+2]) == "{"
    body = line.split("{")[1].split("}")[0].strip()
    declarations = [l.strip() for l in body.split(";")][:-1]
    fields = [d.split(" ")[1].strip() for d in declarations]
    return name, fields

