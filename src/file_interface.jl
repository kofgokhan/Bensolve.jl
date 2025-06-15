function _vlp_init(filename, vlp, opt)
    result = vlp_init(filename, vlp, opt)
    return result, vlp, opt
end

function solve(filename)
    vlp, opt = _get_default_vlp(), _get_default_opt()
    GC.@preserve vlp opt begin
        mktempdir("."; prefix="jl_Bensolve_") do tmp
            ret, vlp, opt = _vlp_init(filename, vlp, opt)
            # Set output location
            set_opt(opt, 5, ["./bensolve", filename, "-s", "--output_filename", tmp * "/"])
            # Build a solution object
            sol = _sol_init(vlp, opt)
            # Prepare to solve
            lp_init(vlp)
            elapsed_time = @elapsed ret = alg(sol, vlp, opt) # need to handle when ret != 0
            if ret > 0
                write_log_file(vlp, sol, opt, Cdouble(elapsed_time), lp_get_num(Csize_t(0)));
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