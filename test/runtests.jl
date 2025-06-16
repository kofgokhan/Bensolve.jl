using Bensolve
using Test

include("examples/example_utils.jl")

@testset "Bensolve.jl using *.vlp files" begin
    for i = 1:6
        solve("examples/ex0$(i).vlp")
    end
end

@testset "Bensolve.jl using functional interface" begin
    for i = 1:6
        filename = "examples/ex0$(i).vlp"

        opt_dir, m, n, q, cone... = _extract_problem_info_from_file(filename)

        if isempty(cone)
            P = _extract_objective_matrix_from_file(filename, q, n)
            B = _extract_constraint_matrix_from_file(filename, m, n)
            a, b = _extract_row_bounds(filename, m)
            l, u = _extract_col_bounds(filename, n)
            molp_solve(P, B, a, b, l, s, opt_dir)
        else
            ctype, n_gen, nzgen = cone
            P = _extract_objective_matrix_from_file(filename, q, n)
            B = _extract_constraint_matrix_from_file(filename, m, n)
            a, b = _extract_row_bounds(filename, m)
            l, u = _extract_col_bounds(filename, n)
            C = _extract_generator_matrix_from_file(filename, q, n_gen)
            c = _extract_duality_vec(filename, q)
            cone_gen = (lowercase(ctype) == "cone" ? 0 : 1)
            vlp_solve(P, B, a, b, l, u, C, c, opt_dir, cone_gen)
        end
    end
end
