module Bensolve

export molp_solve, vlp_solve

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

function vlp_solve(
    P::AbstractMatrix{<:Real},
    B::AbstractMatrix{<:Real},
    a::AbstractVector{<:Real},
    b::AbstractVector{<:Real},
    l::AbstractVector{<:Real},
    s::AbstractVector{<:Real}, 
    generator_matrix::AbstractMatrix{<:Real}, 
    duality_vector::AbstractVector{<:Real}, 
    opt_dir::Int = 1, 
    cone_gen::cone_gen_type = DEFAULT,
)
    # Getting the dimensions.
    m, n, q = size(B)..., size(P, 1)
    nz, nzobj = sum(B .!= 0), sum(P.!= 0)

    # Building A_ext for vlptype
    A_ext = [B zeros(m, q); -P LinearAlgebra.I(q)]
    l2d, l2d_ref = _matrix2list2d(A_ext)

    # Building the rows and cols for vlptype
    rows, row_refs = _vecs2boundlist(a, b)
    cols, col_refs = _vecs2boundlist(l, s)

    n_gen = size(generator_matrix, 2)
    
    gen = eltype(generator_matrix)[]
    if cone_gen in [CONE, DUALCONE]
        gen = collect(Iterators.flatten(generator_matrix'))
    else
        gen = collect(Iterators.flatten(LinearAlgebra.I(n_gen)'))
    end
    gen_ptr = pointer(gen)

    c = duality_vector
    c_ptr = pointer(c)

    # Get a default opttype
    opt = _get_default_opt()
    
    GC.@preserve l2d l2d_ref rows row_refs cols col_refs gen_ptr c_ptr begin
        # Build the problem
        vlp = vlptype(
            pointer_from_objref(l2d),
            pointer_from_objref(rows), 
            pointer_from_objref(cols), 
            opt_dir,
            cone_gen,
            gen_ptr,
            c_ptr,
            nz,
            nzobj,
            n, m, q, n_gen
        )
        # Solve and store the solution temporarily
        mktempdir("."; prefix="jl_Bensolve_") do tmp
            vlp_ref = Ref(vlp)
            # Set output location
            set_opt(opt, 5, ["./bensolve", "", "-s", "--output_filename", tmp * "/"])
            # Build a solution object
            sol = _sol_init(vlp_ref, opt)
            # Prepare to solve
            lp_init(vlp_ref)
            elapsed_time = @elapsed ret = alg(sol, vlp_ref, opt) # need to handle when ret != 0
            dump(sol)
            if ret > 0
                write_log_file(vlp_ref, sol, opt, Cdouble(elapsed_time), lp_get_num(Csize_t(0)));
		        display_info(opt, Cdouble(elapsed_time), lp_get_num(Csize_t(0)));
            end
            status = sol[].status
            @info status
            lp_free(Csize_t(0))
            upper_img = Dict{_SolutionIndex, _Solution}()
            lower_img = Dict{_SolutionIndex, _Solution}()
            if status == VLP_OPTIMAL
                adj = _load_adj_inc_info(tmp * "/_adj_p.sol")
                inc = _load_adj_inc_info(tmp * "/_inc_p.sol")
                open(tmp * "/_img_p.sol") do f_img
                    open(tmp * "/_pre_img_p.sol") do f_pre_img
                        for (i, (line_img, line_pre_img)) in enumerate(zip(eachline(f_img), eachline(f_pre_img)))
                            id = i - 1
                            isvertex, y... = parse.(Float64, split(line_img))
                            x = parse.(Float64, split(line_pre_img))
                            upper_img[_SolutionIndex(id)] = _Solution(id, y, x, adj[id], isvertex == 1)
                        end
                    end
                end
                open(tmp * "/_img_d.sol") do f_img
                    open(tmp * "/_pre_img_d.sol") do f_pre_img
                        for (i, (line_img, line_pre_img)) in enumerate(zip(eachline(f_img), eachline(f_pre_img)))
                            id = i - 1
                            isvertex, y... = parse.(Float64, split(line_img))
                            x = parse.(Float64, split(line_pre_img))
                            lower_img[_SolutionIndex(id)] = _Solution(id, y, x, inc[id], isvertex == 1)
                        end
                    end
                end
            end
            return status, upper_img, lower_img, elapsed_time
        end
    end
end

function molp_solve(
    P::AbstractMatrix{<:Real},
    B::AbstractMatrix{<:Real},
    a::AbstractVector{<:Real},
    b::AbstractVector{<:Real},
    l::AbstractVector{<:Real},
    s::AbstractVector{<:Real}, 
    opt_dir::Int = 1;
)
    # Getting the dimensions.
    m, n, q = size(B)..., size(P, 1)
    nz, nzobj = sum(B .!= 0), sum(P.!= 0)

    # Building A_ext for vlptype
    A_ext = [B zeros(m, q); -P LinearAlgebra.I(q)]
    l2d, l2d_ref = _matrix2list2d(A_ext)

    # Building the rows and cols for vlptype
    rows, row_refs = _vecs2boundlist(a, b)
    cols, col_refs = _vecs2boundlist(l, s)

    # Get a default opttype
    opt = _get_default_opt()

    status = VLP_NOSTATUS
    
    GC.@preserve l2d l2d_ref rows row_refs cols col_refs begin
        # Build the problem
        vlp = vlptype(
            pointer_from_objref(l2d),
            pointer_from_objref(rows), 
            pointer_from_objref(cols), 
            opt_dir,
            DEFAULT,
            C_NULL,
            C_NULL,
            nz,
            nzobj,
            n, m, q, 0
        )
        # Solve and store the solution temporarily
        mktempdir("."; prefix="jl_Bensolve_") do tmp
            vlp_ref = Ref(vlp)
            # Set output location
            set_opt(opt, 5, ["./bensolve", "", "-s", "--output_filename", tmp * "/"])
            # Build a solution object
            sol = _sol_init(vlp_ref, opt)
            # Prepare to solve
            lp_init(vlp_ref)
            elapsed_time = @elapsed ret = alg(sol, vlp_ref, opt) # need to handle when ret != 0
            dump(sol)
            if ret > 0
                write_log_file(vlp_ref, sol, opt, Cdouble(elapsed_time), lp_get_num(Csize_t(0)));
		        display_info(opt, Cdouble(elapsed_time), lp_get_num(Csize_t(0)));
            end
            status = sol[].status
            @info status
            lp_free(Csize_t(0))
            upper_img = Dict{_SolutionIndex, _Solution}()
            lower_img = Dict{_SolutionIndex, _Solution}()
            if status == VLP_OPTIMAL
                adj = _load_adj_inc_info(tmp * "/_adj_p.sol")
                inc = _load_adj_inc_info(tmp * "/_inc_p.sol")
                open(tmp * "/_img_p.sol") do f_img
                    open(tmp * "/_pre_img_p.sol") do f_pre_img
                        for (i, (line_img, line_pre_img)) in enumerate(zip(eachline(f_img), eachline(f_pre_img)))
                            id = i - 1
                            isvertex, y... = parse.(Float64, split(line_img))
                            x = parse.(Float64, split(line_pre_img))
                            upper_img[_SolutionIndex(id)] = _Solution(id, y, x, adj[id], isvertex == 1)
                        end
                    end
                end
                open(tmp * "/_img_d.sol") do f_img
                    open(tmp * "/_pre_img_d.sol") do f_pre_img
                        for (i, (line_img, line_pre_img)) in enumerate(zip(eachline(f_img), eachline(f_pre_img)))
                            id = i - 1
                            isvertex, y... = parse.(Float64, split(line_img))
                            x = parse.(Float64, split(line_pre_img))
                            lower_img[_SolutionIndex(id)] = _Solution(id, y, x, inc[id], isvertex == 1)
                        end
                    end
                end
            end
            return status, upper_img, lower_img, elapsed_time
        end
    end
end

function _load_adj_inc_info(filename)
    dict = Dict{Int, Vector{Int}}()
    open(filename) do f
        for (i, l) in enumerate(eachline(f))
            id = i - 1
            neighbors = parse.(Int, split(l))
            dict[id] = neighbors
        end
    end
    return dict
end

function solve(
    P::AbstractMatrix{<:Real},
    B::AbstractMatrix{<:Real},
    a::AbstractVector{<:Real},
    b::AbstractVector{<:Real},
    l::AbstractVector{<:Real},
    s::AbstractVector{<:Real};
    opt_dir::Int = 1,
    cone_gen::cone_gen_type = DEFAULT,
    generator_matrix::AbstractMatrix{<:Real}, 
    duality_vector::AbstractVector{<:Real}, 
)
    # Getting the dimensions.
    m, n, q = size(B)..., size(P, 1)
    nz, nzobj = sum(B .!= 0), sum(P.!= 0)

    # Building A_ext for vlptype
    A_ext = [B zeros(m, q); -P LinearAlgebra.I(q)]
    l2d, A_ext_refs = _matrix2list2d(A_ext)

    # Building the rows and cols for vlptype
    rows, row_refs = _vecs2boundlist(a, b)
    cols, col_refs = _vecs2boundlist(l, s)

    n_gen = 0
    if generator_matrix !== nothing
        n_gen = size(generator_matrix, 2)
    end

    gen_ptr = C_NULL
    gen = zeros(q * n_gen)
    if generator_matrix !== nothing
        gen .= collect(Iterators.flatten(generator_matrix'))
        gen_ptr = pointer(gen)
    end

    c_ptr = C_NULL
    c = zeros(q)
    if duality_vector !== nothing
        c .= duality_vector
        c_ptr = pointer(c)
    end

    opt = _get_default_opt()

    # GC.@preserve l2d idx1 idx2 data rows row_refs cols col_refs gen c opt begin
    GC.@preserve l2d A_ext_refs rows row_refs cols col_refs gen c opt begin
        vlp = vlptype(
            pointer_from_objref(l2d),
            pointer_from_objref(rows), 
            pointer_from_objref(cols), 
            opt_dir,
            cone_gen,
            gen_ptr,
            c_ptr,
            nz,
            nzobj,
            n, m, q, n_gen
        )
        mktempdir("."; prefix="jl_Bensolve_") do tmp
            vlp_ref = Ref(vlp)
            set_opt(opt, 4, ["./bensolve", "tempfile", "--output_filename", tmp * "/"])
            sol = _sol_init(vlp_ref, opt)
            lp_init(vlp_ref)
            ret = alg(sol, vlp_ref, opt)
            res = readdlm(tmp * "/_img_p.sol")
            vertices = res[res[:, 1] .== 1, 2:end]
            directions = res[res[:, 1] .== 0, 2:end]
            return vertices, directions
        end
    end
end

function _create_molp(
    P::AbstractMatrix{<:Real},
    B::AbstractMatrix{<:Real},
    a::AbstractVector{<:Real},
    b::AbstractVector{<:Real},
    l::AbstractVector{<:Real},
    s::AbstractVector{<:Real};
    opt_dir::Int = 1, 
)
    # Getting the dimensions.
    m, n, q = size(B)..., size(P, 1)
    nz, nzobj = sum(B .!= 0), sum(P.!= 0)

    # Building A_ext for vlptype
    A_ext = [B zeros(m, q); -P LinearAlgebra.I(q)]
    l2d, l2d_ref = _matrix2list2d(A_ext)

    # Building the rows and cols for vlptype
    rows, row_refs = _vecs2boundlist(a, b)
    cols, col_refs = _vecs2boundlist(l, s)

    # Get a default opttype
    opt = _get_default_opt()
    
    GC.@preserve l2d l2d_ref rows row_refs cols col_refs begin
        # Build the problem
        vlp = vlptype(
            pointer_from_objref(l2d),
            pointer_from_objref(rows), 
            pointer_from_objref(cols), 
            opt_dir,
            DEFAULT,
            C_NULL,
            C_NULL,
            nz,
            nzobj,
            n, m, q, 0
        )
    end

    return vlp, (l2d, l2d_ref, rows, row_refs, cols, col_refs)
end

function _create_vlp(
    P::AbstractMatrix{<:Real},
    B::AbstractMatrix{<:Real},
    a::AbstractVector{<:Real},
    b::AbstractVector{<:Real},
    l::AbstractVector{<:Real},
    s::AbstractVector{<:Real};
    opt_dir::Int = 1,
    cone_gen::cone_gen_type = DEFAULT,
    generator_matrix::AbstractMatrix{<:Real}, 
    duality_vector::AbstractVector{<:Real}, 
)
    # Getting the dimensions.
    m, n, q = size(B)..., size(P, 1)
    nz, nzobj = sum(B .!= 0), sum(P.!= 0)

    # Building A_ext for vlptype
    A_ext = [B zeros(m, q); -P LinearAlgebra.I(q)]
    l2d, l2d_ref = _matrix2list2d(A_ext)

    # Building the rows and cols for vlptype
    rows, row_refs = _vecs2boundlist(a, b)
    cols, col_refs = _vecs2boundlist(l, s)

    # Get a default opttype
    opt = _get_default_opt()
    
    GC.@preserve l2d l2d_ref rows row_refs cols col_refs begin
        # Build the problem
        vlp = vlptype(
            pointer_from_objref(l2d),
            pointer_from_objref(rows), 
            pointer_from_objref(cols), 
            opt_dir,
            DEFAULT,
            C_NULL,
            C_NULL,
            nz,
            nzobj,
            n, m, q, 0
        )
    end
end

end