function _get_vlp(
    P::AbstractMatrix, 
    B::AbstractMatrix, 
    a::AbstractVector, 
    b::AbstractVector, 
    l::AbstractVector, 
    s::AbstractVector,  
    opt_dir::Int = 1, 
    cone_gen::cone_gen_type = DEFAULT; 
    Y_or_Z::AbstractMatrix = Float64[], 
)
    m, n, q = size(B)..., size(P, 1)
    n_gen = size(Y_or_Z, 2)
    A_ext = _matrix2list2d([B zeros(m, q); P -I(q)])
    nz, nzobj = nnz(sparse(B)), nnz(sparse(P))
    gen = cone_gen == DEFAULT ? C_NULL : Y_or_Z
    rows = _vecs2rows(a, b)
    cols = _vecs2cols(l, s)

    vlp = Ref(vlptype(
        Ptr{list2d}(A_ext),           # A_ext
        Ptr{boundlist}(rows),         # rows
        Ptr{boundlist}(cols),         # cols
        opt_dir,                      # optdir
        cone_gen_type,                # cone_gen
        Ptr{Cdouble}(gen),            # gen
        Ptr{Cdouble}(zeros(q)),       # c
        Clong(nz),                    # nz
        Clong(nzobj),                 # nzobj
        lp_idx(n),                    # n
        lp_idx(m),                    # m
        lp_idx(q),                    # q
        lp_idx(n_gen)                 # n_gen
    ))
    return vlp
end

function _matrix2list2d(A::AbstractMatrix{<:Real})
    I, J, V = findnz(sparse(A))

    idx1 = lp_idx.(I)
    idx2 = lp_idx.(J)
    data = Cdouble.(V)

    l2d = list2d(
        Csize_t(length(data)),
        pointer(idx1),
        pointer(idx2),
        pointer(data)
    )

    return l2d, (idx1, idx2, data)
end

function _vecs2boundlist(lb::AbstractVector{<:Real}, ub::AbstractVector{<:Real})
    @assert length(lb) == length(ub)
    dim = length(lb)

    idx = lp_idx.(1:dim)
    lb_ = zeros(Cdouble, dim)
    ub_ = zeros(Cdouble, dim)
    typ_ = Cchar.(fill('d', dim))
    for i in 1:dim
        if isinf(lb[i]) && isinf(ub[i])
            typ_[i] = 'f'
        elseif isinf(lb[i])
            typ_[i] = 'u'
            ub_[i] = ub[i]
        elseif isinf(ub[i])
            typ_[i] = 'l'
            lb_[i] = lb[i]
        elseif lb[i] == ub[i]
            typ_[i] = 's'
            lb_[i] = lb[i]
            ub_[i] = ub[i]
        else
            typ_[i] = 'd'
            lb_[i] = lb[i]
            ub_[i] = ub[i]
        end
    end

    bounds = boundlist(
        lp_idx(dim),
        pointer(idx),
        pointer(lb_),
        pointer(ub_),
        pointer(typ_)
    )

    return bounds, (idx, lb_, ub_, typ_)
end

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