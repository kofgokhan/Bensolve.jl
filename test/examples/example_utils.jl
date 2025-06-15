using Bensolve
using SparseArrays

function _extract_problem_info_from_file(filename)
    open(filename) do f
        for line in eachline(f)
            desc, _, data... = split(strip(line))
            if desc == "p"
                sense, vlp... = data
                opt_dir = lowercase(sense) == "min" ? 1 : -1
                m, n, nz, q, nzobj, cone... = something.(tryparse.(Int, vlp), vlp)
                if !isempty(cone)
                    cone_type, n_gen, nz_gen = cone
                    return opt_dir, m, n, nz, q, nzobj, cone_type, n_gen, nz_gen
                end
                return opt_dir, m, n, nz, q, nzobj
            end
        end
    end
end

function _extract_matrix_from_file(filename, designator, nrows, ncols)
    matrix = open(filename) do f
        I = Int[]
        J = Int[]
        V = Float64[]
        for line in eachline(f)
            desc, data... = split(line)
            if desc == designator
                i, j, v = parse.(Float64, data)
                if j > 0
                    push!(I, i)
                    push!(J, j)
                    push!(V, v)
                end
            end
        end
        sparse(I, J, V, nrows, ncols)
    end
    return matrix
end

function _bound_from_desc(typ, bound)
    if typ == "l"
        return (only(bound), Inf)
    elseif typ == "u"
        return (-Inf, only(bound))
    elseif typ == "f"
        return (-Inf, Inf)
    elseif typ == "d"
        return bound
    else
        return (only(bound), only(bound))
    end
end

function _extract_bounds(filename, designator, s)
    l = fill(-Inf, s)
    u = fill(Inf, s)
    open(filename) do f
        for line in eachline(f)
            desc, data... = split(strip(line))
            if desc == designator
                index, typ, bound... = data
                i = parse(Int, index)
                l[i], u[i] = _bound_from_desc(typ, tryparse.(Float64, bound))
            end
        end
    end
    return l, u
end

function _extract_duality_vec(filename, s)
    c = fill(1., s)
    open(filename) do f
        for line in eachline(f)
            desc, data... = split(strip(line))
            if desc == "k"
                i, j, v = parse.([Int, Int, Float64], data)
                if j == 0
                    c[i] = v
                end
            end
        end
    end
    return c
end

_extract_objective_matrix_from_file(filename, nrows, ncols) = _extract_matrix_from_file(filename, "o", nrows, ncols)
_extract_constraint_matrix_from_file(filename, nrows, ncols) = _extract_matrix_from_file(filename, "a", nrows, ncols)
_extract_generator_matrix_from_file(filename, nrows, ncols) = _extract_matrix_from_file(filename, "k", nrows, ncols)

_extract_row_bounds(filename, s) = _extract_bounds(filename, "i", s)
_extract_col_bounds(filename, s) = _extract_bounds(filename, "j", s)

fname = "test/examples/ex05.vlp"

opt_dir, m, n, nz, q, nzobj, ctype, n_gen, nzgen = _extract_problem_info_from_file(fname)

P = _extract_objective_matrix_from_file(fname, q, n)
B = _extract_constraint_matrix_from_file(fname, m, n)
C = _extract_generator_matrix_from_file(fname, q, n_gen)
a, b = _extract_row_bounds(fname, m)
l, s = _extract_col_bounds(fname, n)
c = _extract_duality_vec(fname, q)

Bensolve.molp_solve(P, B, a, b, l, s)
status, upper_img, _ = Bensolve.vlp_solve(P, B, a, b, l, s, C, c, 1, Bensolve.CONE)
upper_img