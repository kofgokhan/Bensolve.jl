function _get_default_opt()
    opt = Ref(
        opttype(
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
            DEFAULT_EPS_BENSON_PHASE2,
        ),
    )
    set_default_opt(opt)
    return opt
end

function _get_default_vlp()
    vlp = Ref(
        vlptype(
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
            lp_idx(0),                     # n_gen
        ),
    )
    return vlp
end

function _get_default_sol()
    return Ref(
        soltype(
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
            Csize_t(0),             # dd_dir
        ),
    )
end

function _sol_init(vlp, opt)
    sol = _get_default_sol()
    sol_init(sol, vlp, opt)
    return sol
end
