function evaluate_individual!(phenotype, individual, optsys)
    newprob = remake_odeprob(individual, optsys)

    saved_values = SavedValues(Float64, Tuple{Float64, Float64})
    cb = SavingCallback(optsys.amem_saver, saved_values; saveat=0.1)

    # Update and solve the ODE problem with the given parameters
    # sol, saved_values = solve_odes(individual, odeprob)
    sol = solve(newprob, optsys.alg, callback=cb)

    if !successful_retcode(sol.retcode)
        print("Unsuccessful")
        return phenotype
    end

    # Reshape the saved values to a 2D array
    saved_array = reinterpret(reshape, Float64, saved_values.saveval)

    # Divide the saved values by the initial/total A to get the normalized values
    # saved_array ./= A_total
    
    # Calculate and return the fitness of the solution
    return calculate_fitness!(phenotype, saved_array, optsys)
end


# Out of place solver, mainly used for testing or solving the fittest individual during the trace output. Takes an individual in the form of a DimVector, creates a SavedValues object to store the results, and returns the solution and the saved values
function solve_odes(individual::DimVector, optsys)
    #- Remake the ODEProblem with the new parameters and species
    newprob = remake_odeprob(individual, optsys)
    #- Create a SavedValues object to store the results
    saved_values = SavedValues(Float64, Tuple{Float64, Float64})
    #- Create a callback to save the results at specified time intervals
    cb = SavingCallback(optsys.amem_saver, saved_values, saveat=0.1)
    #- Solve the ODE problem
    sol = solve(newprob, optsys.alg, callback=cb)

    # Reshape the saved values to a 2D array
    saved_array = reinterpret(reshape, Float64, saved_values.saveval)

    # Divide the saved values by the initial/total A to get the normalized values
    # saved_array ./= A_total

    return sol, saved_array
end

# Remakes the ODEProblem with the new parameters and species. Uses the parameter and species setter functions defined in `OptimizationReactionSystem`.
function remake_odeprob(individual::DimVector, optsys)
    # Divides the individual into parameters and species correctly
    u0, p = get_u0_p(individual)
    # Need to copy the original MTKParameters when remaking in order to avoid issues with mutability
    newprob = remake(optsys.oprob, p = copy(optsys.oprob.p), u0 = copy(optsys.oprob.u0))
    optsys.parameter_setter(newprob, parent(p))
    optsys.species_setter(newprob, parent(u0))
    return newprob
end






