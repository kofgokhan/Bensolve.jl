import MathOptInterface as MOI

@enum(_VariableBound, _LOWER, _UPPER, _DOUBLE, _FREE, _FIXED,)

_bounds(s::MOI.GreaterThan{Float64}) = (s.lower, Inf)
_bounds(s::MOI.LessThan{Float64}) = (-Inf, s.upper)
_bounds(s::MOI.EqualTo{Float64}) = (s.value, s.value)
_bounds(s::MOI.Interval{Float64}) = (s.lower, s.upper)

struct _VariableInfo
    index::MOI.VariableIndex
    column::Int
    bound::_VariableBound
    name::String
end

function _VariableInfo(index::MOI.VariableIndex, column::Int)
    return _VariableInfo(index, column, _FREE, "")
end

mutable struct _ConstraintInfo
    row::Int
    set::MOI.AbstractSet
    name::String
    _ConstraintInfo(set) = new(0, set, "")
    _ConstraintInfo(row, set) = new(row, set, "")
end

struct _SolutionIndex
    value::Int
end

Base.convert(::Type{Int}, x::_SolutionIndex) = x.value
Base.convert(::Type{_SolutionIndex}, i::Int) = _SolutionIndex(i)
Base.:(==)(i::_SolutionIndex, j::_SolutionIndex) = i.value == j.value
Base.hash(idx::_SolutionIndex, h::UInt) = hash(idx.value, h)

struct _Solution
    id::_SolutionIndex
    x::Vector{Float64}
    y::Vector{Float64}
    adj::Vector{Int}
    isvertex::Bool
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    map::MOI.IndexMap
    ctype::cone_gen_type
    upper_img::Dict{_SolutionIndex,_Solution}
    lower_img::Dict{_SolutionIndex,_Solution}
    options::Dict{String,String}
    status::sol_status_type
    solve_time::Float64

    Optimizer() = new(
        MOI.IndexMap(),
        DEFAULT,
        Dict{_SolutionIndex,_Solution}(),
        Dict{_SolutionIndex,_Solution}(),
        Dict{String,String}(),
        VLP_NOSTATUS,
        0.0,
    )
end

function MOI.empty!(optimizer::Optimizer)
    optimizer.map = MOI.IndexMap()
    optimizer.ctype = DEFAULT
    optimizer.upper_img = Dict{_SolutionIndex,_Solution}()
    optimizer.lower_img = Dict{_SolutionIndex,_Solution}()
    optimizer.options = Dict{String,String}()
    optimizer.status = VLP_NOSTATUS
    optimizer.solve_time = 0.0
    return
end

function MOI.is_empty(optimizer::Optimizer)
    return isempty(optimizer.map) &&
           optimizer.ctype == DEFAULT &&
           isempty(optimizer.upper_img) &&
           isempty(optimizer.lower_img) &&
           isempty(optimizer.options) &&
           optimizer.status == VLP_NOSTATUS &&
           optimizer.solve_time == 0.0
end

function MOI.optimize!(dest::Optimizer, src::MOI.ModelLike)
    MOI.empty!(dest)
    sense = MOI.get(src, MOI.ObjectiveSense())
    variables = MOI.get(src, MOI.ListOfVariableIndices())
    n = 0
    for x in variables
        n += 1
        dest.map[x] = MOI.VariableIndex(n)
    end
    l, u = fill(-Inf, n), fill(Inf, n)
    for S in [
        MOI.LessThan{Float64},
        MOI.EqualTo{Float64},
        MOI.GreaterThan{Float64},
        MOI.Interval{Float64},
    ]
        constraints = MOI.get(src, MOI.ListOfConstraintIndices{MOI.VariableIndex,S}())
        for ci in constraints
            f = MOI.get(src, MOI.ConstraintFunction(), ci)
            s = MOI.get(src, MOI.ConstraintSet(), ci)
            l[dest.map[f].value], u[dest.map[f].value] = _bounds(s)
        end
    end

    objs = MOI.Utilities.scalarize(
        MOI.get(src, MOI.ObjectiveFunction{MOI.VectorAffineFunction{Float64}}()),
    )
    q = length(objs)
    objective_matrix = spzeros(Float64, q, n)
    obj_constants = zeros(q)
    for (i, obj) in enumerate(objs)
        obj_constants[i] = obj.constant
        for term in obj.terms
            objective_matrix[i, dest.map[term.variable].value] = term.coefficient
        end
    end

    m = 0
    for S in [
        MOI.LessThan{Float64},
        MOI.EqualTo{Float64},
        MOI.GreaterThan{Float64},
        MOI.Interval{Float64},
    ]
        constraints =
            MOI.get(src, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},S}())
        for ci in constraints
            m += 1
            dest.map[ci] = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S}(m)
        end
    end
    constraint_matrix = spzeros(Float64, m, n)
    a, b = fill(-Inf, m), fill(Inf, m)
    for S in [MOI.LessThan{Float64}, MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}]
        constraints =
            MOI.get(src, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},S}())
        for ci in constraints
            f = MOI.get(src, MOI.ConstraintFunction(), ci)
            for term in f.terms
                constraint_matrix[dest.map[ci].value, dest.map[term.variable].value] =
                    term.coefficient
            end
            s = MOI.get(src, MOI.ConstraintSet(), ci)
            a[dest.map[ci].value], b[dest.map[ci].value] = _bounds(s)
        end
    end
    status, upper_img, lower_img, solve_time = molp_solve(
        objective_matrix,
        constraint_matrix,
        a,
        b,
        l,
        u,
        sense == MOI.MIN_SENSE ? 1 : -1,
    )
    dest.status = status
    dest.upper_img = upper_img
    dest.lower_img = lower_img
    dest.solve_time = solve_time

    return dest.map, false
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Bensolve"

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{<:Union{MOI.LessThan,MOI.GreaterThan,MOI.EqualTo}},
)
    return true
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.VectorAffineFunction{Float64}},
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarAffineFunction{Float64}},
    ::Type{<:Union{MOI.EqualTo{Float64},MOI.LessThan{Float64},MOI.GreaterThan{Float64}}},
)
    return true
end

struct OrderingCone <: MOI.AbstractModelAttribute end

MOI.supports(::Optimizer, ::OrderingCone) = true

function MOI.set(model::Optimizer, ::OrderingCone, C::AbstractMatrix{<:Real})
    model.generator_matrix = C
    return
end

function MOI.get(model::Optimizer, ::OrderingCone)
    return model.generator_matrix
end

struct DualityVector <: MOI.AbstractModelAttribute end

MOI.supports(::Optimizer, ::DualityVector) = true

function MOI.set(model::Optimizer, ::DualityVector, c::AbstractVector{<:Real})
    model.duality_vector = c
    return
end

struct ConeType <: MOI.AbstractModelAttribute end

MOI.supports(::Optimizer, ::ConeType) = true

function MOI.set(model::Optimizer, ::ConeType, ctype::cone_gen_type)
    model.ctype = ctype
    return
end

function MOI.get(model::Optimizer, ::ConeType)
    return model.ctype
end

function MOI.get(model::Optimizer, ::MOI.DualStatus)
    if model.status == VLP_OPTIMAL
        return MOI.FEASIBLE_POINT
    else
        return MOI.NO_SOLUTION
    end
end

function MOI.get(model::Optimizer, ::MOI.PrimalStatus)::MOI.ResultStatusCode
    if model.status == VLP_OPTIMAL
        return MOI.FEASIBLE_POINT
    elseif model.status == VLP_INFEASIBLE
        return MOI.INFEASIBLE_POINT
    else
        return MOI.NO_SOLUTION
    end
end

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    if model.status == VLP_UNBOUNDED
        return "VLP is totally unbounded, there is no solution"
    elseif model.status == VLP_NOVERTEX
        return "upper image of VLP has no vertex (this case is not covered by this version)"
    elseif model.status == VLP_INFEASIBLE
        return "VLP is infeasible"
    elseif model.status == VLP_UNBOUNDED
        return "VLP is not bounded"
    elseif model.status == VLP_OPTIMAL
        return "VLP is solved"
    else
        return "Unknown Status"
    end
end

MOI.get(model::Optimizer, ::MOI.ResultCount) = length(model.upper_img)

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)::MOI.TerminationStatusCode
    if model.status == VLP_OPTIMAL
        return MOI.OPTIMAL
    elseif model.status == VLP_INFEASIBLE
        return MOI.INFEASIBLE
    elseif model.status == VLP_UNBOUNDED
        return MOI.DUAL_INFEASIBLE
    elseif model.status == VLP_NOVERTEX
        return MOI.INFEASIBLE_OR_UNBOUNDED
    elseif model.status in (VLP_INPUTERROR, VLP_UNEXPECTED_STATUS)
        return MOI.OTHER_ERROR
    end
end

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    return model.upper_img[_SolutionIndex(attr.result_index - 1)].y
end

function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    return model.solve_time
end

function MOI.get(model::Optimizer, attr::MOI.VariablePrimal, x::MOI.VariableIndex)
    model.upper_img[_SolutionIndex(attr.result_index - 1)].x[model.map[x].value]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{<:MOI.AbstractScalarFunction},
)
    return model.lower_img[_SolutionIndex(attr.result_index - 1)].y[model.map[ci].value]
end

function MOI.get(model::Optimizer, attr::MOI.DualObjectiveValue)
    return model.lower_img[_SolutionIndex(attr.result_index - 1)].y
end
