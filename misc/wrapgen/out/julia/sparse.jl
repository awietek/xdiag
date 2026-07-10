# SPDX-License-Identifier: Apache-2.0
# STATIC hand-written special -- copied verbatim by generate.sh.
#
# Julia-owned sparse matrices, filled in place by the C++ two-phase kernels.
# int64 / int32 index type is chosen by the function name (matching the C++
# library); the coefficient type (Float64/ComplexF64) is picked at runtime from
# isreal. i0 (index base) is a keyword, as in the C++ library. The structs
# expose their raw arrays, so callers who want a SparseArrays.SparseMatrixCSC
# can build one themselves (the library does not depend on SparseArrays).

export coo_matrix,
    coo_matrix_32, csr_matrix, csr_matrix_32, csc_matrix, csc_matrix_32, to_dense

_sparse_coeff(ops, block) = (isreal(ops) && isreal(block)) ? Float64 : ComplexF64

function _coo_matrix(ops::OpSum, block_in::Block, ::Type{Ti}; i0 = 1) where {Ti<:Integer}
    Tc = _sparse_coeff(ops, block_in)
    herm = ishermitian(ops, block_in)
    nnz = fun_coo_nnz(ops.cxx_object, block_in.cxx_object)
    nrows = Ti(fun_sparse_dim_out(ops.cxx_object, block_in.cxx_object))
    ncols = Ti(size(block_in))
    row, col, data = Vector{Ti}(undef, nnz), Vector{Ti}(undef, nnz), Vector{Tc}(undef, nnz)
    GC.@preserve row col data begin
        fun_coo_fill(
            ops.cxx_object,
            block_in.cxx_object,
            Int64(nnz),
            pointer(row),
            pointer(col),
            pointer(data),
            Ti(i0),
        )
    end
    return COOMatrix{Ti,Tc}(nrows, ncols, row, col, data, Ti(i0), herm)
end
coo_matrix(ops::OpSum, block::Block; i0 = 1) = _coo_matrix(ops, block, Int64; i0)
coo_matrix_32(ops::OpSum, block::Block; i0 = 1) = _coo_matrix(ops, block, Int32; i0)

function _csr_matrix(ops::OpSum, block_in::Block, ::Type{Ti}; i0 = 1) where {Ti<:Integer}
    Tc = _sparse_coeff(ops, block_in)
    herm = ishermitian(ops, block_in)
    counts = fun_csr_nnz(ops.cxx_object, block_in.cxx_object, false)   # per-row
    nrows, ncols, nnz = Ti(length(counts)), Ti(size(block_in)), sum(counts)
    rowptr, col, data =
        Vector{Ti}(undef, nrows + 1), Vector{Ti}(undef, nnz), Vector{Tc}(undef, nnz)
    GC.@preserve rowptr col data begin
        fun_csr_fill(
            ops.cxx_object,
            block_in.cxx_object,
            counts,
            pointer(rowptr),
            pointer(col),
            pointer(data),
            Ti(i0),
            false,
        )
    end
    return CSRMatrix{Ti,Tc}(nrows, ncols, rowptr, col, data, Ti(i0), herm)
end
csr_matrix(ops::OpSum, block::Block; i0 = 1) = _csr_matrix(ops, block, Int64; i0)
csr_matrix_32(ops::OpSum, block::Block; i0 = 1) = _csr_matrix(ops, block, Int32; i0)

function _csc_matrix(ops::OpSum, block_in::Block, ::Type{Ti}; i0 = 1) where {Ti<:Integer}
    Tc = _sparse_coeff(ops, block_in)
    herm = ishermitian(ops, block_in)
    counts = fun_csr_nnz(ops.cxx_object, block_in.cxx_object, true)    # per-column (transpose)
    ncols, nnz = Ti(length(counts)), sum(counts)
    nrows = Ti(fun_sparse_dim_out(ops.cxx_object, block_in.cxx_object))
    colptr, row, data =
        Vector{Ti}(undef, ncols + 1), Vector{Ti}(undef, nnz), Vector{Tc}(undef, nnz)
    GC.@preserve colptr row data begin
        fun_csr_fill(
            ops.cxx_object,
            block_in.cxx_object,
            counts,
            pointer(colptr),
            pointer(row),
            pointer(data),
            Ti(i0),
            true,
        )
    end
    return CSCMatrix{Ti,Tc}(nrows, ncols, colptr, row, data, Ti(i0), herm)
end
csc_matrix(ops::OpSum, block::Block; i0 = 1) = _csc_matrix(ops, block, Int64; i0)
csc_matrix_32(ops::OpSum, block::Block; i0 = 1) = _csc_matrix(ops, block, Int32; i0)

# matvec / matmat y = A*x. A real CSRMatrix applied to a complex x promotes to a
# complex result; a complex CSRMatrix always yields a complex result (a real x is
# promoted). The C++ apply chooses the (possibly MKL) kernel by argument type.
function _csr_apply_eltype(A::CSRMatrix, x)
    if eltype(A.data) <: Real
        return eltype(x) <: Real ? Float64 : ComplexF64
    else
        return ComplexF64
    end
end

apply(A::CSRMatrix, x::AbstractVector) =
    apply(A, x, Vector{_csr_apply_eltype(A, x)}(undef, Int(A.nrows)))
apply(A::CSRMatrix, x::AbstractMatrix) =
    apply(A, x, Matrix{_csr_apply_eltype(A, x)}(undef, Int(A.nrows), size(x, 2)))

function apply(A::CSRMatrix, x::AbstractVector, y::AbstractVector)
    xin = eltype(x) === eltype(y) ? x : convert(Vector{eltype(y)}, x)
    GC.@preserve A xin y fun_csr_apply_vec(
        csr_view(A),
        pointer(xin),
        pointer(y),
        length(xin),
        length(y),
    )
    return y
end

function apply(A::CSRMatrix, x::AbstractMatrix, y::AbstractMatrix)
    xin = eltype(x) === eltype(y) ? x : convert(Matrix{eltype(y)}, x)
    GC.@preserve A xin y fun_csr_apply_mat(
        csr_view(A),
        pointer(xin),
        pointer(y),
        size(xin, 1),
        size(xin, 2),
        size(y, 1),
    )
    return y
end

# dense: all three via the C++ to_dense (CSR through its view; COO/CSC fused).
to_dense(A::CSRMatrix) = GC.@preserve A to_julia(fun_csr_to_dense(csr_view(A)))
to_dense(A::COOMatrix) = GC.@preserve A to_julia(
    fun_coo_to_dense(
        A.nrows,
        A.ncols,
        pointer(A.row),
        pointer(A.col),
        pointer(A.data),
        length(A.data),
        A.i0,
    ),
)
to_dense(A::CSCMatrix) = GC.@preserve A to_julia(
    fun_csc_to_dense(
        A.nrows,
        A.ncols,
        pointer(A.colptr),
        pointer(A.row),
        pointer(A.data),
        length(A.data),
        A.i0,
    ),
)
