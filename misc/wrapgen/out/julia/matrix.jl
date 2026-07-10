# SPDX-License-Identifier: Apache-2.0
# STATIC hand-written special -- copied verbatim by generate.sh.
#
# Dense matrix representation of an OpSum on a block, allocated Julia-side and
# filled by the C++ kernel. Real inputs give a Matrix{Float64}, else complex.

export matrix

function matrix(ops::OpSum, block_in::Block)
    nrows = fun_sparse_dim_out(ops.cxx_object, block_in.cxx_object)   # size(out)
    ncols = dim(block_in)
    if isreal(ops) && isreal(block_in)
        m = Matrix{Float64}(undef, nrows, ncols)
        GC.@preserve m fun_matrix(ops.cxx_object, block_in.cxx_object, pointer(m))
        return m
    else
        m = Matrix{ComplexF64}(undef, nrows, ncols)
        GC.@preserve m fun_matrixC(ops.cxx_object, block_in.cxx_object, pointer(m))
        return m
    end
end

# Op / Monomial: a one-term OpSum, so reuse the zero-copy OpSum path.
matrix(op::Op, block::Block) = matrix(OpSum(op), block)
matrix(mono::Monomial, block::Block) = matrix(OpSum(mono), block)
