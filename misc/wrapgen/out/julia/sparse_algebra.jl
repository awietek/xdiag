# SPDX-License-Identifier: Apache-2.0
# STATIC hand-written special -- copied verbatim by generate.sh.
#
# The Julia-owned sparse matrix types and the C++ view helper. These are a
# foundation module (included before the generated submodules) so the wrapped
# CSRMatrix algorithms can annotate their arguments `::CSRMatrix`. Builders and
# to_dense live in sparse.jl; the algorithms are generated.

struct COOMatrix{Ti<:Integer,Tc<:Number}
    nrows::Ti
    ncols::Ti
    row::Vector{Ti}
    col::Vector{Ti}
    data::Vector{Tc}
    i0::Ti
    ishermitian::Bool
end

struct CSRMatrix{Ti<:Integer,Tc<:Number}
    nrows::Ti
    ncols::Ti
    rowptr::Vector{Ti}
    col::Vector{Ti}
    data::Vector{Tc}
    i0::Ti
    ishermitian::Bool   # copied into the C++ view; drives the apply() matvec mode
end

struct CSCMatrix{Ti<:Integer,Tc<:Number}
    nrows::Ti
    ncols::Ti
    colptr::Vector{Ti}
    row::Vector{Ti}
    data::Vector{Tc}
    i0::Ti
    ishermitian::Bool
end

export COOMatrix, CSRMatrix, CSCMatrix

# Non-owning C++ CSRMatrix view over this matrix's Julia-owned arrays. The
# caller must GC.@preserve the CSRMatrix for the duration of the C++ call (the
# generated algorithm forwarders do this).
csr_view(A::CSRMatrix) = fun_csr_view(
    A.nrows,
    A.ncols,
    pointer(A.rowptr),
    pointer(A.col),
    pointer(A.data),
    length(A.data),
    A.i0,
    A.ishermitian,
)
