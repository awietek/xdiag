class_template_instantiations = \
    {"COOMatrix": {"idx_t": ["int32_t", "int64_t"],
                   "coeff_t": ["double", "complex"]},
     "CSRMatrix": {"idx_t": ["int32_t", "int64_t"],
                   "coeff_t": ["double", "complex"]},
     "CSCMatrix": {"idx_t": ["int32_t", "int64_t"],
                   "coeff_t": ["double", "complex"]}}

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
    
ignored = ["operator<<", "operator++", "operator[]", "fill", "what",
           "irreps", # in blocks: would return internal RepresenationSet
           "basis", # in blocks: would return std::shared_ptr
           "table", # in FileToml: would return the raw toml table
           "begin", "end",
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
                 "Error"]
blocks = ["Fermion", "Boson", "Spinhalf", "Electron", "tJ"]
