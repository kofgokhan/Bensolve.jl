using Revise
using JuMP
using SparseArrays

import MathOptInterface as MOI

using Bensolve

include("examples/example_utils.jl")

filename = "test/examples/ex01.vlp"

function test_problem(filename)
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
        ctype, C, c = cone
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

test_problem("test/examples/ex01.vlp")

filename = "test/examples/ex05.vlp"
opt_dir, m, n, q, cone... = _extract_problem_info_from_file(filename)
P = _extract_objective_matrix_from_file(filename, q, n)
B = _extract_constraint_matrix_from_file(filename, m, n)
a, b = _extract_row_bounds(filename, m)
l, u = _extract_col_bounds(filename, n)
ctype, n_gen, nzgen = cone
C = _extract_generator_matrix_from_file(filename, q, n_gen)
c = _extract_duality_vec(filename, n_gen)

cone_gen = (lowercase(ctype) == "cone" ? 1 : 2)
Bensolve.vlp_solve(P, B, a, b, l, u, C, c, 1, cone_gen)