using Bensolve
using Test

include("examples/example_utils.jl")

@testset "Bensolve.jl using *.vlp files" begin
    for i = 1:6
        0 == solve("examples/ex0$(i).vlp")
    end
end

@testset "Bensolve.jl using functional interface" begin
    for i = 1:6
        filename = "test/examples/ex0$(i).vlp"

        opt_dir, m, n, nz, q, nzobj, cone... = _extract_problem_info_from_file(fname)

        if isempty(cone)
            P = _extract_objective_matrix_from_file(fname, q, n)
            B = _extract_constraint_matrix_from_file(fname, m, n)
            a, b = _extract_row_bounds(fname, m)
            l, s = _extract_col_bounds(fname, n)
            molp_solve(P, B, a, b, l, s, opt_dir)
        else
            ctype, C, c = cone
            P = _extract_objective_matrix_from_file(fname, q, n)
            B = _extract_constraint_matrix_from_file(fname, m, n)
            a, b = _extract_row_bounds(fname, m)
            l, s = _extract_col_bounds(fname, n)
            C = _extract_generator_matrix_from_file(fname, q, n_gen)
            c = _extract_duality_vec(fname, q)
            vlp_solve(P, B, a, b, l, s, opt_dir, ctype, C, c)
        end
    end
end
