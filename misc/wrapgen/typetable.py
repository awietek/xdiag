# SPDX-License-Identifier: Apache-2.0
#
# The type table: the single place that knows, per supported type, its Julia
# name, its C++ spelling, and how to convert a value across the boundary.
#
# Every function here takes a type STRING as it appears in the model, e.g.
# "int64_t", "Op const&", "State&", "double*". decay() drops the qualifier to
# get the table key; the qualifier itself decides what can be wrapped:
#   value / const& / (const& return)  -> wrappable
#   mutable T& / T&& / T*             -> NOT wrappable (skip -> hand-write)

PRIMITIVE = {"bool": "Bool", "int64_t": "Int64", "double": "Float64", "complex": "ComplexF64"}
STRING = {"std::string": "String"}
ARMA = {                                        # key -> (julia, cpp)
    "arma::vec": ("Vector{Float64}", "arma::vec"),
    "arma::cx_vec": ("Vector{ComplexF64}", "arma::cx_vec"),
    "arma::mat": ("Matrix{Float64}", "arma::mat"),
    "arma::cx_mat": ("Matrix{ComplexF64}", "arma::cx_mat"),
    "arma::ivec": ("Vector{Int64}", "arma::Col<int64_t>"),
    "arma::imat": ("Matrix{Int64}", "arma::Mat<int64_t>"),
}
BLOCKS = ["Fermion", "Boson", "Spinhalf", "Electron", "tJ"]

CLASSES = set()                                 # filled from the model by emit.py
RESULT_STRUCTS = set()                           # records with public `fields`


def register_classes(names):
    CLASSES.update(names)


def register_result_structs(names):
    RESULT_STRUCTS.update(names)


def decay(typestr):
    for q in ("const&", "&&", "&", "*"):
        typestr = typestr.replace(q, "")
    return typestr.strip()


def ref_kind(typestr):
    if typestr.endswith("*"):
        return "pointer"
    if typestr.endswith("&&"):
        return "rvalue"
    if typestr.endswith("&"):
        return "const-ref" if "const" in typestr else "mutable-ref"
    return "value"


def category(typestr):
    key = decay(typestr)
    if key == "void":
        return "void"
    if key in PRIMITIVE:
        return "primitive"
    if key in STRING:
        return "string"
    if key in ARMA:
        return "arma"
    if key == "Block":
        return "block"
    if key.startswith("tuple<"):
        return "tuple"
    if key == "CSRMatrix":          # Julia-owned sparse matrix, viewed by C++
        return "sparse"
    if key in CLASSES:
        return "wrapped"
    return "unsupported"


def tuple_elements(typestr):
    inner = decay(typestr)[len("tuple<"):-1]
    return [e.strip() for e in inner.split(",")]


def supported(typestr, as_return=False):
    cat = category(typestr)
    if cat == "unsupported":
        return False
    if cat == "tuple":                          # returned as a Julia tuple only
        return as_return and all(supported(e, as_return=True) for e in tuple_elements(typestr))
    if cat == "sparse":                         # CSRMatrix only as a (viewed) argument
        return not as_return
    if cat == "block" and as_return:            # can't build a Block variant
        return False
    kind = ref_kind(typestr)
    if kind == "pointer":
        return False
    if not as_return and kind == "rvalue":
        return False
    if not as_return and kind == "mutable-ref" and cat != "wrapped":
        return False        # value-type out-param can't propagate; wrapped T& is in-place
    return True


# --- rendering --------------------------------------------------------------
def julia_type(typestr):
    cat, key = category(typestr), decay(typestr)
    if cat == "primitive":
        return PRIMITIVE[key]
    if cat == "string":
        return STRING[key]
    if cat == "arma":
        return ARMA[key][0]
    if cat == "block":
        return "Block"
    if cat == "wrapped":
        return key
    if cat == "sparse":
        return "CSRMatrix"
    if cat == "tuple":
        return "Tuple{" + ", ".join(julia_type(e) for e in tuple_elements(typestr)) + "}"
    return "Nothing"


def cpp_param(typestr, name, block_sub=None, sparse_sub=None):
    """C++ lambda parameter: primitives by value, everything else by const&.
    sparse_sub=(idx,coeff) instantiates a CSRMatrix<idx,coeff> argument."""
    cat, key = category(typestr), decay(typestr)
    if cat == "sparse" and sparse_sub:
        base = f"{key}<{sparse_sub[0]}, {sparse_sub[1]}>"
    elif cat == "block" and block_sub:
        base = block_sub
    elif cat == "arma":
        base = ARMA[key][1]
    else:
        base = key
    if cat == "primitive":
        return f"{base} {name}"
    if ref_kind(typestr) == "mutable-ref":      # in-place: pass the shared object by ref
        return f"{base} &{name}"
    return f"{base} const& {name}"


def unwrap(typestr, expr):
    """Julia value -> argument for the CxxWrap call."""
    cat = category(typestr)
    if cat == "arma":
        return f"to_armadillo({expr})"
    if cat == "sparse":                 # build a non-owning C++ view of the arrays
        return f"csr_view({expr})"
    if cat in ("wrapped", "block"):
        return f"{expr}.cxx_object"
    return expr


def wrap(typestr, expr):
    """CxxWrap result -> Julia value."""
    cat = category(typestr)
    if cat == "arma":
        return f"to_julia({expr})"
    if cat == "string":                         # std::string -> native Julia String
        return f"String({expr})"
    if cat == "wrapped":
        key = decay(typestr)
        if key in RESULT_STRUCTS:               # native struct built from accessors
            return f"convert({key}, {expr})"
        return f"{key}({expr})"
    if cat == "tuple":                          # (r = call; (wrap(r[1]), wrap(r[2]), ...))
        parts = ", ".join(wrap(e, f"r[{i + 1}]") for i, e in enumerate(tuple_elements(typestr)))
        return f"(r = {expr}; ({parts}))"
    return expr


# index-shift helpers (broadcast form works for scalars and arma vectors alike)
def shift_down(expr):
    return f"({expr} .- 1)"


def shift_up(expr):
    return f"({expr} .+ 1)"
