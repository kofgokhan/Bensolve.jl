module Bensolve

export molp_solve, vlp_solve, solve

import libbensolve_jll
using SparseArrays
using LinearAlgebra
using DelimitedFiles

function __init__()
    global libbensolve = libbensolve_jll.libbensolve
    return
end

include("gen/libbensolve.jl")
include("common.jl")
include("file_interface.jl")
include("functional_interface.jl")
include("MOI_wrapper/MOI_wrapper.jl")

function _load_adj_inc_info(filename)
    dict = Dict{Int,Vector{Int}}()
    open(filename) do f
        for (i, l) in enumerate(eachline(f))
            id = i - 1
            neighbors = parse.(Int, split(l))
            dict[id] = neighbors
        end
    end
    return dict
end

end
