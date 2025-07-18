# Out of place solver, mainly used for testing or solving the fittest individual during the trace output. Takes an individual in the form of a DimVector, creates a SavedValues object to store the results, and returns the solution and the saved values
function solve_odes(individual::DimVector, optsys::OptimizationReactionSystem)
    newprob = remake_odeprob(individual, optsys.oprob, optsys.parameter_setter, optsys.species_setter)
    saved_values = SavedValues(Float64, Float64)
    cb = SavingCallback(optsys.amem_saver, saved_values, saveat=0.1)
    sol = solve(newprob, Rodas5P(chunk_size=optsys.autodiff_chunk_size), callback=cb)

    return sol, saved_values
end

# In place solver that is called by `evaluate_individual!` during the optimization process. The SavedValues is created in `evaluate_individual!` and mutated during the solve here.
function solve_odes!(individual::DimVector, integrator, optsys::OptimizationReactionSystem)
    newprob = remake_odeprob(individual, optsys.oprob, optsys.parameter_setter, optsys.species_setter)
    sol = solve(newprob, Rodas5P(chunk_size=optsys.autodiff_chunk_size), callback=cb)
    return sol
end





# Remakes the ODEProblem with the new parameters and species. Uses the parameter and species setter functions defined in `OptimizationReactionSystem`.
function remake_odeprob(individual::DimVector, odeprob, parameter_setter, species_setter)
    # Divides the individual into parameters and species correctly
    u0, p = get_u0_p(individual)
    # Need to copy the original MTKParameters when remaking in order to avoid issues with mutability
    newprob = remake(odeprob, p = copy(odeprob.p), u0 = copy(odeprob.u0))
    parameter_setter(newprob, parent(p))
    species_setter(newprob, parent(u0))
    return newprob
end







# DEPRECATED #
# For analysis of datasets 
function get_u0_params(individual_row::DataFrameRow)
    u0 = NamedTuple(individual_row[Between(:L, :A)]) |> pairs |> collect 
    params = NamedTuple(individual_row[Between(:kfᴸᴬ, :DF)]) |> pairs |> collect
    return u0, params
end

function solve_odes(individual_row::DataFrameRow, odeprob)
    u0, p = get_u0_params(individual_row)
    newprob = remake(odeprob, u0 = u0, p = p)
    return solve(newprob, Rodas4P(chunk_size = length(odeprob.u0)), saveat = 0.1, save_on = true)
end