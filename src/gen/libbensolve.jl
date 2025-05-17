using CEnum

const __jmp_buf = NTuple{8, Clong}

struct __sigset_t
    __val::NTuple{16, Culong}
end

struct __jmp_buf_tag
    __jmpbuf::__jmp_buf
    __mask_was_saved::Cint
    __saved_mask::__sigset_t
end

const jmp_buf = NTuple{1, __jmp_buf_tag}

const lp_idx = Cint

@cenum _pre_img_type::UInt32 begin
    PRE_IMG_OFF = 0
    PRE_IMG_ON = 1
end

const pre_img_type = _pre_img_type

@cenum _format_type::UInt32 begin
    FORMAT_SHORT = 0
    FORMAT_LONG = 1
    FORMAT_AUTO = 2
end

const format_type = _format_type

@cenum _lp_method_type::UInt32 begin
    PRIMAL_SIMPLEX = 0
    DUAL_SIMPLEX = 1
    DUAL_PRIMAL_SIMPLEX = 2
    LP_METHOD_AUTO = 3
end

const lp_method_type = _lp_method_type

@cenum _alg_type::UInt32 begin
    PRIMAL_BENSON = 0
    DUAL_BENSON = 1
end

const alg_type = _alg_type

struct opttype
    bounded::Cint
    plot::Cint
    filename::NTuple{256, Cchar}
    solution::pre_img_type
    format::format_type
    lp_method_phase0::lp_method_type
    lp_method_phase1::lp_method_type
    lp_method_phase2::lp_method_type
    message_level::Cint
    lp_message_level::Cint
    alg_phase1::alg_type
    alg_phase2::alg_type
    eps_phase0::Cdouble
    eps_phase1::Cdouble
    eps_benson_phase1::Cdouble
    eps_benson_phase2::Cdouble
end

@cenum _cone_out_type::UInt32 begin
    CONE_OUT_OFF = 0
    CONE_OUT_ON = 1
end

const cone_out_type = _cone_out_type

@cenum _swap_type::UInt32 begin
    SWAP = 0
    NO_SWAP = 1
end

const swap_type = _swap_type

function cone_vertenum(prim, n_prim, dual, n_dual, prim_in, n_prim_in, dim, opt, output, swap)
    ccall((:cone_vertenum, libbensolve), Cint, (Ptr{Ptr{Cdouble}}, Ptr{lp_idx}, Ptr{Ptr{Cdouble}}, Ptr{lp_idx}, Ptr{Cdouble}, Csize_t, Csize_t, Ptr{opttype}, cone_out_type, swap_type), prim, n_prim, dual, n_dual, prim_in, n_prim_in, dim, opt, output, swap)
end

@cenum _sol_status_type::UInt32 begin
    VLP_NOSTATUS = 0
    VLP_INFEASIBLE = 1
    VLP_UNBOUNDED = 2
    VLP_NOVERTEX = 3
    VLP_OPTIMAL = 4
    VLP_INPUTERROR = 5
    VLP_UNEXPECTED_STATUS = 6
end

const sol_status_type = _sol_status_type

@cenum _c_dir_type::UInt32 begin
    C_DIR_POS = 0
    C_DIR_NEG = 1
end

const c_dir_type = _c_dir_type

struct soltype
    m::lp_idx
    n::lp_idx
    q::lp_idx
    o::lp_idx
    p::lp_idx
    r::lp_idx
    h::lp_idx
    eta::Ptr{Cdouble}
    Y::Ptr{Cdouble}
    Z::Ptr{Cdouble}
    c::Ptr{Cdouble}
    R::Ptr{Cdouble}
    H::Ptr{Cdouble}
    status::sol_status_type
    c_dir::c_dir_type
    pp::Csize_t
    dd::Csize_t
    pp_dir::Csize_t
    dd_dir::Csize_t
end

mutable struct list2d
    size::Csize_t
    idx1::Ptr{lp_idx}
    idx2::Ptr{lp_idx}
    data::Ptr{Cdouble}
end

mutable struct boundlist
    size::lp_idx
    idx::Ptr{lp_idx}
    lb::Ptr{Cdouble}
    ub::Ptr{Cdouble}
    type::Ptr{Cchar}
end

@cenum _cone_gen_type::UInt32 begin
    CONE = 0
    DUALCONE = 1
    DEFAULT = 2
end

const cone_gen_type = _cone_gen_type

mutable struct vlptype
    A_ext::Ptr{list2d}
    rows::Ptr{boundlist}
    cols::Ptr{boundlist}
    optdir::Cint
    cone_gen::cone_gen_type
    gen::Ptr{Cdouble}
    c::Ptr{Cdouble}
    nz::Clong
    nzobj::Clong
    n::lp_idx
    m::lp_idx
    q::lp_idx
    n_gen::lp_idx
end

function alg(sol, vlp, opt)
    ccall((:alg, libbensolve), Cint, (Ptr{soltype}, Ptr{vlptype}, Ptr{opttype}), sol, vlp, opt)
end

function phase0(sol, vlp, opt)
    ccall((:phase0, libbensolve), Cvoid, (Ptr{soltype}, Ptr{vlptype}, Ptr{opttype}), sol, vlp, opt)
end

function phase1_primal(sol, vlp, opt)
    ccall((:phase1_primal, libbensolve), Cvoid, (Ptr{soltype}, Ptr{vlptype}, Ptr{opttype}), sol, vlp, opt)
end

function phase2_primal(sol, vlp, opt)
    ccall((:phase2_primal, libbensolve), Cvoid, (Ptr{soltype}, Ptr{vlptype}, Ptr{opttype}), sol, vlp, opt)
end

function phase1_dual(sol, vlp, opt)
    ccall((:phase1_dual, libbensolve), Cvoid, (Ptr{soltype}, Ptr{vlptype}, Ptr{opttype}), sol, vlp, opt)
end

function phase2_dual(sol, vlp, opt)
    ccall((:phase2_dual, libbensolve), Cvoid, (Ptr{soltype}, Ptr{vlptype}, Ptr{opttype}), sol, vlp, opt)
end

function phase2_init(sol, vlp)
    ccall((:phase2_init, libbensolve), Cvoid, (Ptr{soltype}, Ptr{vlptype}), sol, vlp)
end

struct list1d
    size::lp_idx
    idx::Ptr{lp_idx}
    data::Ptr{Cdouble}
end

function string_fprint(filename, string)
    ccall((:string_fprint, libbensolve), Cvoid, (Ptr{Cchar}, Ptr{Cchar}), filename, string)
end

function matrix_fprint(mat_arr, m, n, tda, filename, format)
    ccall((:matrix_fprint, libbensolve), Cvoid, (Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cchar}, Ptr{Cchar}), mat_arr, m, n, tda, filename, format)
end

function matrix_print(mat_arr, m, n, format)
    ccall((:matrix_print, libbensolve), Cvoid, (Ptr{Cdouble}, Cint, Cint, Ptr{Cchar}), mat_arr, m, n, format)
end

function string_to_int(str, error_msg)
    ccall((:string_to_int, libbensolve), Cint, (Ptr{Cchar}, Ptr{Cchar}), str, error_msg)
end

function string_to_positive_int(str, error_msg)
    ccall((:string_to_positive_int, libbensolve), Cint, (Ptr{Cchar}, Ptr{Cchar}), str, error_msg)
end

function string_to_positive_double(str, error_msg)
    ccall((:string_to_positive_double, libbensolve), Cdouble, (Ptr{Cchar}, Ptr{Cchar}), str, error_msg)
end

function orthogonal_vector(mat_arr, dim, cidx)
    ccall((:orthogonal_vector, libbensolve), Cvoid, (Ptr{Cdouble}, Cint, Cint), mat_arr, dim, cidx)
end

function is_equal(size, vec1, vec2, tol)
    ccall((:is_equal, libbensolve), Cint, (lp_idx, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble), size, vec1, vec2, tol)
end

function is_zero(size, vec, tol)
    ccall((:is_zero, libbensolve), Cint, (lp_idx, Ptr{Cdouble}, Cdouble), size, vec, tol)
end

function list1d_alloc(size)
    ccall((:list1d_alloc, libbensolve), Ptr{list1d}, (lp_idx,), size)
end

function list1d_calloc(size)
    ccall((:list1d_calloc, libbensolve), Ptr{list1d}, (lp_idx,), size)
end

function list1d_init_idx(list, firstidx)
    ccall((:list1d_init_idx, libbensolve), Cvoid, (Ptr{list1d}, lp_idx), list, firstidx)
end

function list1d_free(list)
    ccall((:list1d_free, libbensolve), Cvoid, (Ptr{list1d},), list)
end

function vector_to_list1d(list, vec_arr, n)
    ccall((:vector_to_list1d, libbensolve), Cvoid, (Ptr{list1d}, Ptr{Cdouble}, Cint), list, vec_arr, n)
end

function list1d_print(list, size)
    ccall((:list1d_print, libbensolve), Cvoid, (Ptr{list1d}, Cint), list, size)
end

function list2d_alloc(size)
    ccall((:list2d_alloc, libbensolve), Ptr{list2d}, (Csize_t,), size)
end

function list2d_calloc(size)
    ccall((:list2d_calloc, libbensolve), Ptr{list2d}, (Csize_t,), size)
end

function list2d_init_idx(list, nrows, ncols)
    ccall((:list2d_init_idx, libbensolve), Cvoid, (Ptr{list2d}, lp_idx, lp_idx), list, nrows, ncols)
end

function list2d_free(list)
    ccall((:list2d_free, libbensolve), Cvoid, (Ptr{list2d},), list)
end

function list2d_print(list)
    ccall((:list2d_print, libbensolve), Cvoid, (Ptr{list2d},), list)
end

function boundlist_alloc(size)
    ccall((:boundlist_alloc, libbensolve), Ptr{boundlist}, (lp_idx,), size)
end

function boundlist_calloc(size, type)
    ccall((:boundlist_calloc, libbensolve), Ptr{boundlist}, (lp_idx, Cchar), size, type)
end

function boundlist_init_idx(list, firstidx)
    ccall((:boundlist_init_idx, libbensolve), Cvoid, (Ptr{boundlist}, lp_idx), list, firstidx)
end

function boundlist_free(list)
    ccall((:boundlist_free, libbensolve), Cvoid, (Ptr{boundlist},), list)
end

function boundlist_print(list)
    ccall((:boundlist_print, libbensolve), Cvoid, (Ptr{boundlist},), list)
end

function lp_write(i)
    ccall((:lp_write, libbensolve), Cvoid, (Csize_t,), i)
end

function lp_write_sol(i)
    ccall((:lp_write_sol, libbensolve), Cdouble, (Csize_t,), i)
end

function lp_init(vlp)
    ccall((:lp_init, libbensolve), Cvoid, (Ptr{vlptype},), vlp)
end

@cenum _phase_type::UInt32 begin
    PHASE0 = 0
    PHASE1_PRIMAL = 1
    PHASE1_DUAL = 2
    PHASE2_PRIMAL = 3
    PHASE2_DUAL = 4
end

const phase_type = _phase_type

function lp_set_options(opt, phase)
    ccall((:lp_set_options, libbensolve), Cvoid, (Ptr{opttype}, phase_type), opt, phase)
end

function lp_copy(dest, src)
    ccall((:lp_copy, libbensolve), Cvoid, (Csize_t, Csize_t), dest, src)
end

function lp_update_extra_coeffs(n_rows, n_cols)
    ccall((:lp_update_extra_coeffs, libbensolve), Cvoid, (lp_idx, lp_idx), n_rows, n_cols)
end

function lp_set_rows(i, rows)
    ccall((:lp_set_rows, libbensolve), Cvoid, (Csize_t, Ptr{boundlist}), i, rows)
end

function lp_set_rows_hom(i, rows)
    ccall((:lp_set_rows_hom, libbensolve), Cvoid, (Csize_t, Ptr{boundlist}), i, rows)
end

function lp_set_cols(i, cols)
    ccall((:lp_set_cols, libbensolve), Cvoid, (Csize_t, Ptr{boundlist}), i, cols)
end

function lp_set_cols_hom(i, cols)
    ccall((:lp_set_cols_hom, libbensolve), Cvoid, (Csize_t, Ptr{boundlist}), i, cols)
end

function lp_set_mat(i, A)
    ccall((:lp_set_mat, libbensolve), Cint, (Csize_t, Ptr{list2d}), i, A)
end

function lp_set_mat_row(i, list, ridx)
    ccall((:lp_set_mat_row, libbensolve), Cvoid, (Csize_t, Ptr{list1d}, lp_idx), i, list, ridx)
end

function lp_clear_obj_coeffs(i)
    ccall((:lp_clear_obj_coeffs, libbensolve), Cvoid, (Csize_t,), i)
end

function lp_set_obj_coeffs(i, list)
    ccall((:lp_set_obj_coeffs, libbensolve), Cvoid, (Csize_t, Ptr{list1d}), i, list)
end

@cenum _lp_status_type::UInt32 begin
    LP_INFEASIBLE = 0
    LP_UNBOUNDED = 1
    LP_UNEXPECTED_STATUS = 2
    LP_UNDEFINED_STATUS = 3
    LP_OPTIMAL = 4
end

const lp_status_type = _lp_status_type

function lp_solve(i)
    ccall((:lp_solve, libbensolve), lp_status_type, (Csize_t,), i)
end

function lp_primal_solution_rows(i, x, firstidx, size, sign)
    ccall((:lp_primal_solution_rows, libbensolve), Cvoid, (Csize_t, Ptr{Cdouble}, lp_idx, lp_idx, Cdouble), i, x, firstidx, size, sign)
end

function lp_primal_solution_cols(i, x, firstidx, size, sign)
    ccall((:lp_primal_solution_cols, libbensolve), Cvoid, (Csize_t, Ptr{Cdouble}, lp_idx, lp_idx, Cdouble), i, x, firstidx, size, sign)
end

function lp_dual_solution_rows(i, u, firstidx, size, sign)
    ccall((:lp_dual_solution_rows, libbensolve), Cvoid, (Csize_t, Ptr{Cdouble}, lp_idx, lp_idx, Cdouble), i, u, firstidx, size, sign)
end

function lp_dual_solution_cols(i, u, firstidx, size, sign)
    ccall((:lp_dual_solution_cols, libbensolve), Cvoid, (Csize_t, Ptr{Cdouble}, lp_idx, lp_idx, Cdouble), i, u, firstidx, size, sign)
end

function lp_obj_val(i)
    ccall((:lp_obj_val, libbensolve), Cdouble, (Csize_t,), i)
end

function lp_get_time(i)
    ccall((:lp_get_time, libbensolve), Cdouble, (Csize_t,), i)
end

function lp_get_num(i)
    ccall((:lp_get_num, libbensolve), Cint, (Csize_t,), i)
end

function lp_free(i)
    ccall((:lp_free, libbensolve), Cvoid, (Csize_t,), i)
end

@cenum _lp_hom_type::UInt32 begin
    HOMOGENEOUS = 0
    INHOMOGENEOUS = 1
end

const lp_hom_type = _lp_hom_type

const btstrg = Csize_t

const vrtx_strg = btstrg

struct poly_list_strct
    cnt::Csize_t
    blcks::Csize_t
    data::Ptr{Csize_t}
end

const poly_list = poly_list_strct

struct polytope_strct
    dim::Csize_t
    dim_primg::Csize_t
    cnt::Csize_t
    blcks::Csize_t
    ip::Ptr{Cdouble}
    data::Ptr{Cdouble}
    data_primg::Ptr{Cdouble}
    adjacence::Ptr{poly_list}
    incidence::Ptr{poly_list}
    ideal::Ptr{vrtx_strg}
    used::Ptr{vrtx_strg}
    sltn::Ptr{vrtx_strg}
    dual::Ptr{polytope_strct}
    v2h::Ptr{Cvoid}
end

const polytope = polytope_strct

struct var"##Ctag#257"
    H::Ptr{Cdouble}
    R::Ptr{Cdouble}
    alph::Ptr{Cdouble}
    queue::poly_list
    gnrtrs::poly_list
    intlsd::Cuint
end
function Base.getproperty(x::Ptr{var"##Ctag#257"}, f::Symbol)
    f === :H && return Ptr{Ptr{Cdouble}}(x + 0)
    f === :R && return Ptr{Ptr{Cdouble}}(x + 8)
    f === :alph && return Ptr{Ptr{Cdouble}}(x + 16)
    f === :queue && return Ptr{poly_list}(x + 24)
    f === :gnrtrs && return Ptr{poly_list}(x + 48)
    f === :intlsd && return (Ptr{Cuint}(x + 72), 0, 1)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#257", f::Symbol)
    r = Ref{var"##Ctag#257"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#257"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#257"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


struct poly_args
    data::NTuple{392, UInt8}
end

function Base.getproperty(x::Ptr{poly_args}, f::Symbol)
    f === :dim && return Ptr{Csize_t}(x + 0)
    f === :dim_primg_prml && return Ptr{Csize_t}(x + 8)
    f === :dim_primg_dl && return Ptr{Csize_t}(x + 16)
    f === :ideal && return (Ptr{Cuint}(x + 24), 0, 1)
    f === :idx && return Ptr{Csize_t}(x + 32)
    f === :val && return Ptr{Ptr{Cdouble}}(x + 40)
    f === :val_primg_prml && return Ptr{Ptr{Cdouble}}(x + 48)
    f === :val_primg_dl && return Ptr{Ptr{Cdouble}}(x + 56)
    f === :eps && return Ptr{Cdouble}(x + 64)
    f === :primal && return Ptr{polytope}(x + 72)
    f === :dual && return Ptr{polytope}(x + 184)
    f === :primalV2dualH && return Ptr{Ptr{Cvoid}}(x + 296)
    f === :dualV2primalH && return Ptr{Ptr{Cvoid}}(x + 304)
    f === :init_data && return Ptr{var"##Ctag#257"}(x + 312)
    return getfield(x, f)
end

function Base.getproperty(x::poly_args, f::Symbol)
    r = Ref{poly_args}(x)
    ptr = Base.unsafe_convert(Ptr{poly_args}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{poly_args}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

struct permutation
    cnt::Csize_t
    data::Ptr{Csize_t}
    inv::Ptr{Csize_t}
end

function poly__initialise_permutation(arg1, arg2)
    ccall((:poly__initialise_permutation, libbensolve), Cvoid, (Ptr{polytope}, Ptr{permutation}), arg1, arg2)
end

function poly__kill_permutation(arg1)
    ccall((:poly__kill_permutation, libbensolve), Cvoid, (Ptr{permutation},), arg1)
end

function poly__vrtx2file(arg1, arg2, arg3, arg4)
    ccall((:poly__vrtx2file, libbensolve), Cvoid, (Ptr{polytope}, Ptr{permutation}, Ptr{Cchar}, Ptr{Cchar}), arg1, arg2, arg3, arg4)
end

function poly__primg2file(arg1, arg2, arg3, arg4)
    ccall((:poly__primg2file, libbensolve), Cvoid, (Ptr{polytope}, Ptr{permutation}, Ptr{Cchar}, Ptr{Cchar}), arg1, arg2, arg3, arg4)
end

function poly__adj2file(arg1, arg2, arg3, arg4)
    ccall((:poly__adj2file, libbensolve), Cvoid, (Ptr{polytope}, Ptr{permutation}, Ptr{Cchar}, Ptr{Cchar}), arg1, arg2, arg3, arg4)
end

function poly__inc2file(arg1, arg2, arg3, arg4, arg5)
    ccall((:poly__inc2file, libbensolve), Cvoid, (Ptr{polytope}, Ptr{permutation}, Ptr{permutation}, Ptr{Cchar}, Ptr{Cchar}), arg1, arg2, arg3, arg4, arg5)
end

function poly__set_default_args(args, dim)
    ccall((:poly__set_default_args, libbensolve), Cvoid, (Ptr{poly_args}, Csize_t), args, dim)
end

function poly__initialise(arg1)
    ccall((:poly__initialise, libbensolve), Cvoid, (Ptr{poly_args},), arg1)
end

function poly__kill(arg1)
    ccall((:poly__kill, libbensolve), Cvoid, (Ptr{poly_args},), arg1)
end

function poly__cut(arg1, arg2, arg3)
    ccall((:poly__cut, libbensolve), Cvoid, (Ptr{polytope}, Csize_t, Ptr{Cdouble}), arg1, arg2, arg3)
end

function poly__poly_initialise(arg1, arg2, arg3, arg4, arg5)
    ccall((:poly__poly_initialise, libbensolve), Cvoid, (Ptr{polytope}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Csize_t}), arg1, arg2, arg3, arg4, arg5)
end

function poly__add_vrtx(arg1)
    ccall((:poly__add_vrtx, libbensolve), Cint, (Ptr{poly_args},), arg1)
end

function poly__get_vrtx(arg1)
    ccall((:poly__get_vrtx, libbensolve), Cint, (Ptr{poly_args},), arg1)
end

function poly__poly_init(arg1)
    ccall((:poly__poly_init, libbensolve), Cvoid, (Ptr{polytope},), arg1)
end

function poly__poly_kill(arg1)
    ccall((:poly__poly_kill, libbensolve), Cvoid, (Ptr{polytope},), arg1)
end

function poly__list_init(arg1)
    ccall((:poly__list_init, libbensolve), Cvoid, (Ptr{poly_list},), arg1)
end

function add_vrtx(arg1)
    ccall((:add_vrtx, libbensolve), Cvoid, (Ptr{polytope},), arg1)
end

function poly__intl_apprx(arg1)
    ccall((:poly__intl_apprx, libbensolve), Cint, (Ptr{poly_args},), arg1)
end

function add_lst_elem(arg1, arg2)
    ccall((:add_lst_elem, libbensolve), Cvoid, (Ptr{poly_list}, Csize_t), arg1, arg2)
end

function edge_test(arg1, arg2, arg3)
    ccall((:edge_test, libbensolve), Cint, (Ptr{polytope}, Csize_t, Csize_t), arg1, arg2, arg3)
end

function poly__update_adjacence(arg1)
    ccall((:poly__update_adjacence, libbensolve), Cvoid, (Ptr{polytope},), arg1)
end

function vrtx_cpy(arg1, arg2, arg3)
    ccall((:vrtx_cpy, libbensolve), Cvoid, (Ptr{polytope}, Csize_t, Csize_t), arg1, arg2, arg3)
end

function poly__swap(arg1, arg2)
    ccall((:poly__swap, libbensolve), Cvoid, (Ptr{poly_args}, Ptr{poly_args}), arg1, arg2)
end

function poly__plot(arg1, arg2)
    ccall((:poly__plot, libbensolve), Cvoid, (Ptr{polytope}, Ptr{Cchar}), arg1, arg2)
end

function poly__polyck(arg1)
    ccall((:poly__polyck, libbensolve), Cvoid, (Ptr{poly_args},), arg1)
end

function bslv__normalise(arg1, arg2, arg3, arg4, arg5)
    ccall((:bslv__normalise, libbensolve), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Csize_t, Csize_t), arg1, arg2, arg3, arg4, arg5)
end

struct csatype
    jump::jmp_buf
    fname::Ptr{Cchar}
    fp::Ptr{Libc.FILE}
    count::Cint
    c::Cint
    field::NTuple{256, Cchar}
    empty::Cint
    nonint::Cint
    msg::NTuple{256, Cchar}
    error::Cint
    warning::Cint
end

function vlp_init(filename, vlp, opt)
    ccall((:vlp_init, libbensolve), Cint, (Ptr{Cchar}, Ptr{vlptype}, Ptr{opttype}), filename, vlp, opt)
end

function set_opt(opt, argc, argv)
    ccall((:set_opt, libbensolve), Cint, (Ptr{opttype}, Cint, Ptr{Ptr{Cchar}}), opt, argc, argv)
end

function write_log_file(vlp, sol, opt, elapsedTime, lp_num)
    ccall((:write_log_file, libbensolve), Cint, (Ptr{vlptype}, Ptr{soltype}, Ptr{opttype}, Cdouble, Cint), vlp, sol, opt, elapsedTime, lp_num)
end

function display_info(opt, elapsedTime, lp_num)
    ccall((:display_info, libbensolve), Cvoid, (Ptr{opttype}, Cdouble, Cint), opt, elapsedTime, lp_num)
end

function vlp_free(vlp)
    ccall((:vlp_free, libbensolve), Cvoid, (Ptr{vlptype},), vlp)
end

function sol_init(sol, vlp, opt)
    ccall((:sol_init, libbensolve), Cint, (Ptr{soltype}, Ptr{vlptype}, Ptr{opttype}), sol, vlp, opt)
end

function sol_free(sol)
    ccall((:sol_free, libbensolve), Cvoid, (Ptr{soltype},), sol)
end

function set_default_opt(opt)
    ccall((:set_default_opt, libbensolve), Cvoid, (Ptr{opttype},), opt)
end

const THISVERSION = "2.1.0"

const WELCOME = "BENSOLVE: VLP Solver, %s\nCopyright (C) 2014-2017 Andreas L%shne and Benjamin Wei%sing\nThis is free software; see the source code for copying conditions.\nThere is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or\nFITNESS FOR A PARTICULAR PURPOSE.\n"

const UMLAUT_OE = "ö"

const UMLAUT_SZ = "ß"

const USAGE = "\nUsage: bensolve file [options]\n\n"

const POLY_TEST = 0

const EPS_C = 1.0e-7

const EPS_POLY = 1.0e-9

const EPS_OUTPUT_CHOP = 1.0e-10

const PRIMAL_PLOT_CUT_SHIFT = 0.5

const DUAL_PLOT_CUT_SHIFT = 0.2

const FORMAT_SHORT_STR = "%10.4g "

const FORMAT_LONG_STR = "%.14g "

const DEFAULT_MESSAGE_LEVEL = 1

const DEFAULT_MESSAGE_LEVEL_STR = "1"

const DEFAULT_LP_MESSAGE_LEVEL = 1

const DEFAULT_LP_MESSAGE_LEVEL_STR = "1"

const DEFAULT_EPS_PHASE0 = 1.0e-8

const DEFAULT_EPS_PHASE1 = 1.0e-8

const DEFAULT_EPS_BENSON_PHASE1 = 1.0e-7

const DEFAULT_EPS_BENSON_PHASE1_STR = "1e-7"

const DEFAULT_EPS_BENSON_PHASE2 = 1.0e-7

const DEFAULT_EPS_BENSON_PHASE2_STR = "1e-7"

const MIN_P_STR = "Upper image of primal problem:\n"

const MIN_D_STR = "Lower image of dual problem:\n"

const MAX_P_STR = "Lower image of primal problem:\n"

const MAX_D_STR = "Upper image of dual problem:\n"

const CONE_P_STR = "Ordering cone:\n"

const CONE_D_STR = "Dual of ordering cone:\n"

const CONE_ENDING_STR = ".cone"

const SOL_ENDING_STR = ".sol"

const MAX_STR_LNGTH = 10

const PRE_IMG_P_STR = "_pre_img_p"

const PRE_IMG_D_STR = "_pre_img_d"

const IMG_P_STR = "_img_p"

const IMG_D_STR = "_img_d"

const ADJ_P_STR = "_adj_p"

const ADJ_D_STR = "_adj_d"

const INC_P_STR = "_inc_p"

const INC_D_STR = "_inc_d"

const ALLOCFCTR = 1

# Skipping MacroDefinition: BTCNT ( CHAR_BIT * sizeof ( btstrg ) )

# const VRTXBLCK = ALLOCFCTR * BTCNT

const LSTBLCK = 1

const POLY_EPS = 1.0e-9

const PROBLEM_DESIGNATOR = "vlp"

const STOP_AT_WARNING = 0