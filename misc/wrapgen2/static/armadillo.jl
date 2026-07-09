# SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0

const _ArmadilloArray = Union{Matrix{Int64},
                              Matrix{Float64},
                              Matrix{ComplexF64},
                              Vector{Int64},
                              Vector{Float64},
                              Vector{ComplexF64}}

function _armadillo_wrapper(mat::Matrix{Int64}, copy::Bool)
    m, n = size(mat)
    return typ_arma_mat_int64_t(pointer(mat), m, n, copy, true)
end

function _armadillo_wrapper(mat::Matrix{Float64}, copy::Bool)
    m, n = size(mat)
    return typ_arma_mat(pointer(mat), m, n, copy, true)
end

function _armadillo_wrapper(mat::Matrix{ComplexF64}, copy::Bool)
    m, n = size(mat)
    return typ_arma_cx_mat(pointer(mat), m, n, copy, true)
end

function _armadillo_wrapper(vec::Vector{Int64}, copy::Bool)
    return typ_arma_vec_int64_t(pointer(vec), length(vec), copy, true)
end

function _armadillo_wrapper(vec::Vector{Float64}, copy::Bool)
    return typ_arma_vec(pointer(vec), length(vec), copy, true)
end

function _armadillo_wrapper(vec::Vector{ComplexF64}, copy::Bool)
    return typ_arma_cx_vec(pointer(vec), length(vec), copy, true)
end

function to_armadillo(array::_ArmadilloArray; copy=true)
    GC.@preserve array begin
        return _armadillo_wrapper(array, copy)
    end
end

function with_armadillo(f::F, array::_ArmadilloArray; copy=true) where F
    GC.@preserve array begin
        return f(_armadillo_wrapper(array, copy))
    end
end

function _unsafe_copy_to_julia!(dest::Vector{T}, src, n::Integer) where T
    GC.@preserve dest src begin
        Base.unsafe_copyto!(pointer(dest), memptr(src).cpp_object, n)
    end
    return dest
end

function _unsafe_copy_to_julia!(dest::Matrix{T}, src, n::Integer) where T
    GC.@preserve dest src begin
        Base.unsafe_copyto!(pointer(dest), memptr(src).cpp_object, n)
    end
    return dest
end

function to_julia(vec::typ_arma_vec)
    m = n_rows(vec)
    return _unsafe_copy_to_julia!(Vector{Float64}(undef, m), vec, m)
end

function to_julia(vec::typ_arma_cx_vec)
    m = n_rows(vec)
    return _unsafe_copy_to_julia!(Vector{ComplexF64}(undef, m), vec, m)
end

function to_julia(mat::typ_arma_mat)
    m = n_rows(mat)
    n = n_cols(mat)
    return _unsafe_copy_to_julia!(Matrix{Float64}(undef, m, n), mat, n_elem(mat))
end

function to_julia(mat::typ_arma_cx_mat)
    m = n_rows(mat)
    n = n_cols(mat)
    return _unsafe_copy_to_julia!(Matrix{ComplexF64}(undef, m, n), mat, n_elem(mat))
end
