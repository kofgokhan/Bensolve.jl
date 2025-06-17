using Test
using Bensolve
using JuMP

include("examples/example_utils.jl")

@testset "Bensolve.jl using *.vlp files" begin
    for i = 1:6
        filename = "examples/ex0$(i).vlp"
        solve(filename)
    end
end

@testset "Bensolve.jl using functional interface" begin
    for i = 1:6
        filename = "examples/ex0$(i).vlp"
        cone_gen = 2
        opt_dir, m, n, q, cone... = _extract_problem_info_from_file(filename)

        if isempty(cone)
            P = _extract_objective_matrix_from_file(filename, q, n)
            B = _extract_constraint_matrix_from_file(filename, m, n)
            a, b = _extract_row_bounds(filename, m)
            l, u = _extract_col_bounds(filename, n)
            molp_solve(P, B, a, b, l, u, opt_dir)
        else
            ctype, n_gen = cone
            P = _extract_objective_matrix_from_file(filename, q, n)
            B = _extract_constraint_matrix_from_file(filename, m, n)
            a, b = _extract_row_bounds(filename, m)
            l, u = _extract_col_bounds(filename, n)
            C = _extract_generator_matrix_from_file(filename, q, n_gen)
            c = _extract_duality_vec(filename, q)
            if lowercase(ctype) == "cone"
                cone_gen = 0
            elseif lowercase(ctype) == "dualcone"
                cone_gen = 1
            end
            vlp_solve(P, B, a, b, l, u, C, c, opt_dir, cone_gen)
        end
    end
end


@testset "Bensolve.jl using MOI interface" begin
    for i = 1:6
        filename = "examples/ex0$(i).vlp"
        opt_dir, m, n, q, cone... = _extract_problem_info_from_file(filename)
        P = _extract_objective_matrix_from_file(filename, q, n)
        B = _extract_constraint_matrix_from_file(filename, m, n)
        a, b = _extract_row_bounds(filename, m)
        l, u = _extract_col_bounds(filename, n)
        model = Model(Bensolve.Optimizer)
        @variable(model, l[i] <= x[i=1:n] <= u[i])
        @constraint(model, a .<= B * x .<= b)
        sense = opt_dir == 1 ? MOI.MIN_SENSE : MOI:MAX_SENSE
        @objective(model, sense, P * x)
        if !isempty(cone)
            ctype, n_gen = cone
            C = _extract_generator_matrix_from_file(filename, q, n_gen)
            c = _extract_duality_vec(filename, q)
            MOI.set(model, Bensolve.OrderingCone, C)
            MOI.set(model, Bensolve.DualityVector, c)
            MOI.set(model, Bensolve.ConeType, if ctype == "CONE"
                Bensolve.CONE
            elseif ctype == "DUALCONE"
                Bensolve.DUALCONE
            else
                Bensolve.DEFAULT
            end)
        end
        optimize!(model)
        return
    end
end