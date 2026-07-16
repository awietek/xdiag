# SPDX-License-Identifier: Apache-2.0
# STATIC hand-written special -- copied verbatim by generate.sh.
#
# low-level operator multiplication routine
export apply


function apply(ops::OpSum,
               bin::Block, vin::Vector{Float64},
               bout::Block, vout::Vector{Float64})
    if (length(vin) != length(bin))
        error("Input vector does not have the same size as input Block")
    elseif (length(vout) != length(bout))
        error("Output vector does not have the same size as output Block")
    end
    
    GC.@preserve vin vout begin
        fun_apply(ops.cxx_object, bin.cxx_object, pointer(vin),
                  bout.cxx_object, pointer(vout))
    end
end
function apply(ops::OpSum,
               bin::Block, vin::Vector{ComplexF64},
               bout::Block, vout::Vector{ComplexF64})
    if (length(vin) != length(bin))
        error("Input vector does not have the same size as input Block")
    elseif (length(vout) != length(bout))
        error("Output vector does not have the same size as output Block")
    end
    GC.@preserve vin vout begin
        fun_applyC(ops.cxx_object, bin.cxx_object, pointer(vin),
                                 bout.cxx_object, pointer(vout))
    end
end
function apply(ops::OpSum,
               bin::Block, min::Matrix{Float64},
               bout::Block, mout::Matrix{Float64})
    if (size(min, 1) != length(bin))
        error("Input matrix does not have the same number of rows as the size of the input Block")
    elseif (size(mout, 1) != length(bout))
        error("Output matrix does not have the same number of rows as the size of the output Block")
    elseif size(min, 2) != size(mout, 2)
        error("Input and output matrices need to have the same number of columns")
    end
    GC.@preserve min mout begin
        fun_apply(ops.cxx_object, bin.cxx_object, pointer(min),
                  bout.cxx_object, pointer(mout), size(min, 2))
    end
end
function apply(ops::OpSum,
               bin::Block, min::Matrix{ComplexF64},
               bout::Block, mout::Matrix{ComplexF64})
    if (size(min, 1) != length(bin))
        error("Input matrix does not have the same number of rows as the size of the input Block")
    elseif (size(mout, 1) != length(bout))
        error("Output matrix does not have the same number of rows as the size of the output Block")
    elseif size(min, 2) != size(mout, 2)
        error("Input and output matrices need to have the same number of columns")
    end
    GC.@preserve min mout begin
        fun_applyC(ops.cxx_object, bin.cxx_object, pointer(min),
                   bout.cxx_object, pointer(mout), size(min, 2))
    end
end

# Op / Monomial: a one-term OpSum, so reuse the zero-copy OpSum path.
apply(op::Op, bin::Block, vin, bout::Block, vout) =
    apply(OpSum(op), bin, vin, bout, vout)
apply(op::Monomial, bin::Block, vin, bout::Block, vout) =
    apply(OpSum(op), bin, vin, bout, vout)
