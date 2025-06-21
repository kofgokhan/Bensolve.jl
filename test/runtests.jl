using Test
using Bensolve
using JuMP
using JSON

include("examples/example_utils.jl")

@testset "Bensolve.jl using *.vlp files" begin
    for i = 1:6
        filename = "examples/ex0$(i).vlp"
        status, upper_img, lower_img, solve_time = solve(filename)
        X = [v.x for v in values(upper_img)]
        Y = [v.y for v in values(upper_img)]
        sol = JSON.parsefile("examples/solutions/ex0$(i).json")
        X_true = [Float64.(sol_i["x"]) for sol_i in sol["upper_image"]]
        Y_true = [Float64.(sol_i["y"]) for sol_i in sol["upper_image"]]
        @test sort(X) == sort(X_true)
        @test sort(Y) == sort(Y_true)
    end
end

@testset "Bensolve.jl using functional interface" begin
    for i = 1:6
        filename = "examples/ex0$(i).vlp"
        cone_gen = 2
        opt_dir, m, n, q, cone... = _extract_problem_info_from_file(filename)
        X = []
        Y = []
        if isempty(cone)
            P = _extract_objective_matrix_from_file(filename, q, n)
            B = _extract_constraint_matrix_from_file(filename, m, n)
            a, b = _extract_row_bounds(filename, m)
            l, u = _extract_col_bounds(filename, n)
            status, upper_img, lower_img, solve_time = molp_solve(P, B, a, b, l, u, opt_dir)
            X = [v.x for v in values(upper_img)]
            Y = [v.y for v in values(upper_img)]
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
            status, upper_img, lower_img, solve_time = vlp_solve(P, B, a, b, l, u, C, c, opt_dir, cone_gen)
            X = [v.x for v in values(upper_img)]
            Y = [v.y for v in values(upper_img)]
        end
        sol = JSON.parsefile("examples/solutions/ex0$(i).json")
        X_true = [Float64.(sol_i["x"]) for sol_i in sol["upper_image"]]
        Y_true = [Float64.(sol_i["y"]) for sol_i in sol["upper_image"]]
        @test sort(X) == sort(X_true)
        @test sort(Y) == sort(Y_true)
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
        sense = opt_dir == 1 ? MOI.MIN_SENSE : MOI.MAX_SENSE
        @objective(model, sense, P * x)
        if !isempty(cone)
            ctype, n_gen = cone
            C = _extract_generator_matrix_from_file(filename, q, n_gen)
            c = _extract_duality_vec(filename, q)
            MOI.set(model, Bensolve.OrderingCone(), C)
            MOI.set(model, Bensolve.DualityVector(), c)
            MOI.set(model, Bensolve.ConeType(), if lowercase(ctype) == "cone"
                0
            elseif lowercase(ctype) == "dualcone"
                1
            else
                2
            end)
        end
        optimize!(model)
        X = [value.(x, result=i) for i = 1:result_count(model)]
        Y = [objective_value(model, result=i) for i = 1:result_count(model)]
        sol = JSON.parsefile("examples/solutions/ex0$(i).json")
        X_true = [Float64.(sol_i["x"]) for sol_i in sol["upper_image"]]
        Y_true = [Float64.(sol_i["y"]) for sol_i in sol["upper_image"]]
        @test sort(X) == sort(X_true)
        @test sort(Y) == sort(Y_true)
    end
end