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

# function solve(
#     P::Matrix{Float64}, 
#     B::Matrix{Float64}, 
#     a::Vector{Float64}, 
#     b::Vector{Float64}, 
#     l::Vector{Float64}, 
#     s::Vector{Float64};
#     opt_dir::Int = 1, 
#     cone_gen::cone_gen_type = DEFAULT, 
#     Y_or_Z::Matrix{Float64} = Matrix{Float64}(undef, 2, 2), 
# )
#     m, n, q = size(B)..., size(P, 1)
#     n_gen = size(Y_or_Z, 2)
#     A_ext, A_ext_refs = _matrix2list2d([B zeros(m, q); P -I(q)])
#     nz, nzobj = nnz(sparse(B)), nnz(sparse(P))
#     gen = cone_gen == DEFAULT ? C_NULL : Y_or_Z
#     rows, row_refs = _vecs2rows(a, b)
#     cols, col_refs = _vecs2cols(l, s)
#     c = zeros(q)
# 
#     vlp = Ref(vlptype(
#         Ref(A_ext),                   # A_ext
#         Ref(rows),                    # rows
#         Ref(cols),                    # cols
#         opt_dir,                      # optdir
#         cone_gen_type,                # cone_gen
#         Ref(gen),                     # gen
#         Ref(c),                       # c
#         Clong(nz),                    # nz
#         Clong(nzobj),                 # nzobj
#         lp_idx(n),                    # n
#         lp_idx(m),                    # m
#         lp_idx(q),                    # q
#         lp_idx(n_gen)                 # n_gen
#     ))
# 
#     # vlp, opt = _get_vlp(P, B, a, b, l, s, opt_dir, cone_gen; Y_or_Z), _get_default_opt()
#     GC.@preserve vlp opt begin
#         mktempdir("."; prefix="jl_Bensolve_") do tmp
#             set_opt(opt, 4, ["./bensolve", filename, "--output_filename", tmp * "/"])
#             ret, vlp, opt = _vlp_init(filename, vlp, opt)
#             sol = _sol_init(vlp, opt)
#             lp_init(vlp)
#             ret = alg(sol, vlp, opt)
#             res = readdlm(tmp * "/_img_p.sol")
#             vertices = res[res[:, 1] .== 1, 2:end]
#             directions = res[res[:, 1] .== 0, 2:end]
#             return vertices, directions
#         end
#     end
# end

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