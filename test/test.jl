using Revise
using JuMP
using SparseArrays

import MathOptInterface as MOI

using Bensolve

include("examples/example_utils.jl")

filename = "test/examples/ex06.vlp"

status, upper_img, lower_img, _ = Bensolve.solve(filename)

function test_problem(filename)
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
    return X, Y
end

X, Y = test_problem("test/examples/ex06.vlp")

opt_dir, m, n, q, cone... = _extract_problem_info_from_file(filename)
P = _extract_objective_matrix_from_file(filename, q, n)
B = _extract_constraint_matrix_from_file(filename, m, n)
a, b = _extract_row_bounds(filename, m)
l, u = _extract_col_bounds(filename, n)
ctype, n_gen = cone
C = _extract_generator_matrix_from_file(filename, q, n_gen)
c = _extract_duality_vec(filename, n_gen)

cone_gen = if lowercase(ctype) == "cone" 0 elseif lowercase(ctype) == "dualcone" 1 else 2 end

status, upper_img, _ = Bensolve.vlp_solve(P, B, a, b, l, u, C, c, opt_dir, cone_gen);
