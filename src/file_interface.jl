function _vlp_init(filename, vlp, opt)
    result = vlp_init(filename, vlp, opt)
    return result, vlp, opt
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