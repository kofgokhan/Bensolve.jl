module Bensolve

export _vlp_init, _get_default_opt, _sol_init, solve, alg

import libbensolve_jll
using SparseArrays
using LinearAlgebra
using DelimitedFiles

function __init__()
    global libbensolve = libbensolve_jll.libbensolve
    return
end

include("gen/libbensolve.jl")

function _get_default_opt()
    opt = Ref(opttype(
        0, 
        0, 
        ntuple(_ -> Cchar(0), 256), 
        PRE_IMG_OFF,
        FORMAT_AUTO,
        PRIMAL_SIMPLEX, 
        LP_METHOD_AUTO, 
        LP_METHOD_AUTO,
        DEFAULT_MESSAGE_LEVEL, 
        DEFAULT_LP_MESSAGE_LEVEL,
        PRIMAL_BENSON, 
        PRIMAL_BENSON,
        DEFAULT_EPS_PHASE0, 
        DEFAULT_EPS_PHASE1, 
        DEFAULT_EPS_BENSON_PHASE1, 
        DEFAULT_EPS_BENSON_PHASE2
    ))
    set_default_opt(opt)
    return opt
end

function _get_default_vlp()
    vlp = Ref(vlptype(
        Ptr{list2d}(C_NULL),          # A_ext
        Ptr{boundlist}(C_NULL),       # rows
        Ptr{boundlist}(C_NULL),       # cols
        0,                            # optdir
        CONE,                         # cone_gen
        Ptr{Cdouble}(C_NULL),         # gen
        Ptr{Cdouble}(C_NULL),         # c
        Clong(0),                     # nz
        Clong(0),                     # nzobj
        lp_idx(0),                    # n
        lp_idx(0),                    # m
        lp_idx(0),                    # q
        lp_idx(0)                     # n_gen
    ))
    return vlp
end

function _vlp_init(filename, vlp, opt)
    result = vlp_init(filename, vlp, opt)
    return result, vlp, opt
end

function _get_default_sol()
    return Ref(soltype(
        0,  # m
        0,  # n
        0,  # q
        0,  # o
        0,  # p
        0,  # r
        0,  # h
        Ptr{Cdouble}(C_NULL),  # eta
        Ptr{Cdouble}(C_NULL),  # Y
        Ptr{Cdouble}(C_NULL),  # Z
        Ptr{Cdouble}(C_NULL),  # c
        Ptr{Cdouble}(C_NULL),  # R
        Ptr{Cdouble}(C_NULL),  # H
        VLP_NOSTATUS,          # status (assuming sol_status_type is an enum or primitive)
        C_DIR_POS,             # c_dir (same assumption)
        Csize_t(0),            # pp
        Csize_t(0),            # dd
        Csize_t(0),            # pp_dir
        Csize_t(0)             # dd_dir
    ))
end

function _sol_init(vlp, opt)
    sol = _get_default_sol()
    sol_init(sol, vlp, opt)
    return sol
end

function solve(filename)
    vlp, opt = _get_default_vlp(), _get_default_opt()
    GC.@preserve vlp opt begin
        mktempdir("."; prefix="jl_Bensolve_") do tmp
            set_opt(opt, 4, ["./bensolve", filename, "--output_filename", tmp * "/"])
            ret, vlp, opt = _vlp_init(filename, vlp, opt)
            sol = _sol_init(vlp, opt)
            lp_init(vlp)
            ret = alg(sol, vlp, opt)
            res = readdlm(tmp * "/_img_p.sol")
            vertices = res[res[:, 1] .== 1, 2:end]
            directions = res[res[:, 1] .== 0, 2:end]
            return vertices, directions
        end
    end
end

function solve(
    P::AbstractMatrix, 
    B::AbstractMatrix, 
    a::AbstractVector, 
    b::AbstractVector, 
    l::AbstractVector, 
    s::AbstractVector; 
    Y::AbstractMatrix, 
    Z::AbstractMatrix, 
)
    
end

function _matrix2list2d(A::AbstractMatrix{<:Real})
    A_sparse = sparse(A)
    I, J, V = findnz(A_sparse)

    idx1 = Vector{lp_idx}(lp_idx.(I .- 1))
    idx2 = Vector{lp_idx}(lp_idx.(J .- 1))
    data = Vector{Cdouble}(Cdouble.(V))

    l2d = list2d(
        Csize_t(length(data)),
        pointer(idx1),
        pointer(idx2),
        pointer(data)
    )

    return l2d
end

function _vector2list1d(v::AbstractVector{<:Real})
    v_sparse = sparse(v)
    I, V = findnz(v_sparse)

    idx = Vector{lp_idx}(lp_idx.(I .- 1))
    data = Vector{Cdouble}(Cdouble.(V))

    l1d = list1d(
        lp_idx(length(V)),
        pointer(idx),
        pointer(data)
    )

    return l1d
end

function _vecs2rows(lb::AbstractVector{<:Real}, ub::AbstractVector{<:Real})
    @assert length(lb) == length(ub)
    m = length(lb)

    idx = Vector{lp_idx}(lp_idx.(0:m-1))
    lb_ = Vector{Cdouble}(Cdouble.(lb))
    ub_ = Vector{Cdouble}(Cdouble.(ub))
    typ = Vector{Cchar}(fill('d', m))

    rows = boundlist(
        lp_idx(m),
        pointer(idx),
        pointer(lb_),
        pointer(ub_),
        pointer(typ)
    )

    return rows
end

function _vecs2cols(lb::AbstractVector{<:Real}, ub::AbstractVector{<:Real})
    @assert length(lb) == length(ub)
    n = length(lb)

    idx = Vector{lp_idx}(lp_idx.(0:n-1))
    lb_ = Vector{Cdouble}(Cdouble.(lb))
    ub_ = Vector{Cdouble}(Cdouble.(ub))
    typ = Vector{Cchar}(fill('d', n))

    cols = boundlist(
        lp_idx(n),
        pointer(idx),
        pointer(lb_),
        pointer(ub_),
        pointer(typ)
    )

    return cols
end

end