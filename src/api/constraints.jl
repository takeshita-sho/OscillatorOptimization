# ConstraintFunction represents a single constraint with its function and bounds
struct ConstraintFunction
    func::Function  # The constraint function to be evaluated
    min::Float64    # Lower bound of the constraint
    max::Float64    # Upper bound of the constraint
end

# Custom show method for ConstraintFunction
function Base.show(io::IO, cf::ConstraintFunction)
    print(io, "ConstraintFunction($(cf.min) ≤ f(x) ≤ $(cf.max))")
end

# ConstraintSet holds a collection of ConstraintFunctions
struct ConstraintSet
    constraints::Vector{ConstraintFunction}
end

# Custom show method for ConstraintSet
function Base.show(io::IO, cs::ConstraintSet)
    print(io, "ConstraintSet with $(length(cs.constraints)) constraint(s)")
end

"""
Creates a constraint function for the Km relationship (Kmᴸᴷ < Kmᴸᴾ).
Returns a ConstraintFunction with bounds (-Inf, 0.0).
"""
function km_constraint()
    func = function(individual::DimVector)
        Kmᴸᴷ = (individual[At(:krᴸᴷ)] + individual[At(:kcatᴸᴷ)])/individual[At(:kfᴸᴷ)]
        Kmᴸᴾ = (individual[At(:krᴸᴾ)] + individual[At(:kcatᴸᴾ)])/individual[At(:kfᴸᴾ)]
        return Kmᴸᴷ - Kmᴸᴾ  # Constraint is satisfied when this is negative
    end
    return ConstraintFunction(func, -Inf, 0.0)
end

"""
Creates a constraint function for the K-P relationship (K < P).
Returns a ConstraintFunction with bounds (-Inf, 0.0).
"""
function kp_constraint()
    func = (individual::DimVector) -> individual[At(:K)] - individual[At(:P)]
    return ConstraintFunction(func, -Inf, 0.0)
end

"""
Creates a null constraint that always returns 0.0.
Used as a default when no constraints are specified.
"""
function null_constraint_func(_)
    return 0.0
end

null_constraint() = ConstraintFunction(null_constraint_func, -Inf, Inf)

is_null_constraint(constraint::ConstraintFunction) = constraint.func == null_constraint_func
is_null_constraint(constraint_set::ConstraintSet) = length(constraint_set.constraints) == 1 && is_null_constraint(constraint_set.constraints[1])

"""
Constructs a WorstFitnessConstraints object from a ConstraintSet.
This object is used by the optimization algorithm to enforce constraints.

Args:
- rx_sys: The reaction system
- fixed_params_dict: Dictionary of fixed parameters
- constraint_set: Set of constraints to apply (default: null constraint)

Returns:
WorstFitnessConstraints object with lower/upper bounds and constraint functions.
"""
function make_WorstFitnessConstraints(rx_sys::ReactionSystem, fixed_params_dict::Dict{Symbol, Float64} = Dict{Symbol, Float64}(), constraint_set::ConstraintSet = ConstraintSet([null_constraint()]))

    unfixed_tunable_symbols = get_unfixed_tunable_symbols(fixed_params_dict, rx_sys)
    tunable_bounds_dictionary = get_tunable_bounds_dictionary(unfixed_tunable_symbols, rx_sys)

    # Each entry in tunable_bounds_dictionary is a pair of bounds for a tunable parameter
    # Extract the lower and upper bounds for each tunable parameter into two vectors
    lb = [bounds[1] for bounds in tunable_bounds_dictionary]
    ub = [bounds[2] for bounds in tunable_bounds_dictionary]
    
    # Create a vector of constraint functions
    c = (ind) -> [cf.func(ind) for cf in constraint_set.constraints]

    # Create a vector of constraint lower bounds
    constraint_min = [cf.min for cf in constraint_set.constraints]

    # Create a vector of constraint upper bounds
    constraint_max = [cf.max for cf in constraint_set.constraints]

    # Return a WorstFitnessConstraints object with the lower bounds, upper bounds, constraint lower bounds, constraint upper bounds, and the constraint functions
    return WorstFitnessConstraints(lb, ub, constraint_min, constraint_max, c)
end