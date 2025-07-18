"""
    insert_fixed_inputs!(chromosome::Vector{Pair{Symbol, Float64}}, fixed_inputs_info)

Insert the fixed inputs into the chromosome vector of `Pair{Symbol, Float64}`s, since fixed inputs are not included in the population or it's individuals during optimization.
"""
function insert_fixed_inputs!(chromosome::Vector{Pair{Symbol, Float64}}, fixed_inputs_info)
    for (name, value_and_idx) in fixed_inputs_info
        fixed_value, column_idx = value_and_idx
        insert!(chromosome, column_idx, name => fixed_value)
    end
end

"""Trace override function"""
function trace!(record::Dict{String, Any}, objfun, state, population::SlicedDimArray, method::QD, options) 
    
    #- Get the indices of the oscillatory individuals by checking if the period is greater than 0.0 and amplitude is greater than 0.01
    objective_values = get_objective_values(state)
    oscillatory_population_idxs = findall(col -> col[2] > 0.0 && col[3] > 0.01, eachcol(objective_values))

    #- Get the max fitness value
    record["max_fitness_value"] = value(state)
    @info "Max Fitness Value: ", record["max_fitness_value"]
    
    #- Get the fittest chromosome/individual
    fittest_chromosome = state.fittestChromosome
    record["fittest_individual"] = copy(fittest_chromosome)

    #- Record the oscillatory population array
    oscillatory_population = deepcopy(parent(population)[:, oscillatory_population_idxs])
    record["population"] = oscillatory_population
    println(typeof(oscillatory_population))

    oscillatory_objective_values = objective_values[:, oscillatory_population_idxs]

    # Assert that none of the oscillatory objective values are zeros
    @assert all(oscillatory_objective_values .> 0.0) "Oscillatory objective values are zero"

    record["objective_values"] = oscillatory_objective_values

    # record["fitvals"] = oscillatory_objective_values[1,:]
    # record["periods"] = oscillatory_objective_values[2,:]
    # record["amplitudes"] = oscillatory_objective_values[3,:]

    record["time_vec"] = copy(method.optsys.t)
    # Keep a reference to the optimisation system so that the pretty-print display can
    # access metadata such as the chosen fitness observables.
    record["optsys"] = method.optsys


    # Get Genes symbols 
    genes_symbols = dims(oscillatory_population, :Genes).val.data
    # Add the fixed inputs columns to the results DataFrame
    fixed_inputs_info = oscillatory_population.metadata.val.fixed_inputs_column_info

    # Compute the amem SavedValues of the fittest OSCILLATORY individual here in the trace function, and save it to the record
    # If there are oscillatory solutions, compute and save the solution of the fittest oscillatory individual 
    if length(oscillatory_population_idxs) > 0
        fittest_oscillatory_chromosome = oscillatory_population[:, argmax(oscillatory_objective_values[1,:])]
        fittest_oscillatory_chromosome_map = Pair.(genes_symbols, collect(fittest_oscillatory_chromosome))
        # Add the fixed inputs to the fittest oscillatory chromosome map
        insert_fixed_inputs!(fittest_oscillatory_chromosome_map, fixed_inputs_info)

        println("Fittest oscillatory chromosome: ", fittest_oscillatory_chromosome_map)

        _, saved_vals = solve_odes(fittest_oscillatory_chromosome, method.optsys)
        record["fittest_oscillatory_solution"] = saved_vals
        record["fittest_oscillatory_individual"] = fittest_oscillatory_chromosome
    else
        fittest_chromosome_map = Pair.(genes_symbols, collect(fittest_chromosome))
        insert_fixed_inputs!(fittest_chromosome_map, fixed_inputs_info)
        println("Fittest chromosome: ", fittest_chromosome_map)

        _, saved_vals = solve_odes(state.fittestChromosome, method.optsys)
        record["fittest_solution"] = saved_vals
    end
end




function Evolutionary.show(io::IO, t::Evolutionary.OptimizationTraceRecord{Float64, O}) where {O <: AbstractOptimizer}

    #- Get the max fitness value
    max_fitness_value = t.metadata["max_fitness_value"]

    #- Structure the data for printing
    data = [
        "ITERATION" t.iteration;
        "MAX FITNESS" round(max_fitness_value;digits =2);
        "NUM OSCILLATORY" length(t.metadata["population"]);
    ]

    #- Define the highlighters
    hl1 = Highlighter((data,i,j) -> (i > 1 && j == 2), Crayon(foreground = :white, bold = :true))
    hl2 = Highlighter((data,i,j) -> ([i,j] == [3,2] && data[i,j] == 0), Crayon(foreground = :white, background = :red, bold = :true))
    # hl2 = Highlighter((data,i,j) -> ([i,j] == [3,2] && data[i,j] > 0), Crayon(foreground = :green, background = :black, bold = :true))    

    #- Print the table using PrettyTables.jl
    pretty_table(io, data, header = ["Trace", "Value"], highlighters = (hl1, hl2))


    fittest_individual_key = "fittest_individual"
    fittest_solution_key = "fittest_solution"

    if haskey(t.metadata, "fittest_oscillatory_individual")    
        #- Get fittest oscillatory individual
        fittest_individual_key = "fittest_oscillatory_individual"
        fittest_solution_key = "fittest_oscillatory_solution"
    end
    
    fittest_individual = t.metadata[fittest_individual_key]
    # _, saved_vals = solve_odes(fittest_individual, t.metadata["optsys"])

    # _, saved_vals = solve_odes(maxfit_params, t.metadata["optsys"])
    # Amem = saved_vals.saveval

    #- Get the SavedValues solution object from the trace metadata
    Amem_saved_array = t.metadata[fittest_solution_key]

    fft_obs_sym, td_obs_sym = t.metadata["optsys"].obs_symbols

    obs_fft = Amem_saved_array[1, :]
    obs_td  = Amem_saved_array[2, :]

    # Get new Amem (with AKL and APLp) from the second row of the solution array
    # Amem = Amem_saved_array[1,:]
    # TrimerYield = Amem_saved_array[2,:]


    #- Find the peaks in the Amem solution
    # indx_max, indx_min = findextrema(Amem; min_prominence=0.01)
    amem_indx_max, amem_indx_min = find_amem_peaks(obs_fft)
    amem_vals_max = @view obs_fft[amem_indx_max]
    amem_vals_min = @view obs_fft[amem_indx_min]

    #- Find the peaks in the TrimerYield solution
    trimer_indx_max, trimer_indx_min = find_amem_peaks(obs_td)
    trimer_vals_max = @view obs_td[trimer_indx_max]
    trimer_vals_min = @view obs_td[trimer_indx_min]

    timepoints = t.metadata["time_vec"]

    #- Plot FFT observable (row 1)
    time_plot = lineplot(timepoints, obs_fft; title=fittest_solution_key, xlabel="Time", ylabel="Amplitude", xlim = (minimum(timepoints), maximum(timepoints)), ylim = (0.0, 1.0), color = :red, width= 100, name = string(fft_obs_sym))
    scatterplot!(time_plot, timepoints[amem_indx_max], amem_vals_max; marker = :xcross, color = :blue)
    scatterplot!(time_plot, timepoints[amem_indx_min], amem_vals_min; marker = :xcross, color = :green)
    
    #- Plot time-domain observable (row 2)
    lineplot!(time_plot, timepoints, obs_td; color = :green, name = string(td_obs_sym))
    scatterplot!(time_plot, timepoints[trimer_indx_max], trimer_vals_max; marker = :xcross, color = :blue)
    scatterplot!(time_plot, timepoints[trimer_indx_min], trimer_vals_min; marker = :xcross, color = :green)


    #- Find the peaks in the FFT plot
    fftData = getFrequencies(obs_fft)

    fft_peakindexes = find_fft_peaks(fftData)
    fft_peakvals = @view fftData[fft_peakindexes]
    fft_plot = lineplot(fftData; xlim= (0, maximum(fft_peakindexes)), color = :blue, width= 100, title="FFT")

    scatterplot!(fft_plot, fft_peakindexes, fft_peakvals; marker = :xcross, color = :red)

    periods = t.metadata["objective_values"][2,:]
    amplitudes = t.metadata["objective_values"][3,:]

    #- Histogram of periods 
    if length(periods) < 2
        panel(time_plot) / panel(fft_plot) |> display
        # panel(time_plot) |> display

    else
        period_histogram = histogram(periods; title="Period Histogram", xlabel="Period", color = :green, width= 30)

        #- Histogram of amplitudes
        amplitude_histogram = histogram(amplitudes; title="Amplitude Histogram", xlabel="Amplitude", color = :yellow, width= 30)
    
        panel(time_plot) / panel(fft_plot) / (panel(period_histogram) * panel(amplitude_histogram)) |> display
        # panel(time_plot) / (panel(period_histogram) * panel(amplitude_histogram)) |> display
    end


    println("max(fft obs): ", maximum(obs_fft))
    println("min(fft obs): ", minimum(obs_fft))
    println(fittest_individual_key, ": ", fittest_individual)

    @info "Number of Amem peaks: $(length(amem_indx_max))\nNumber of Amem troughs: $(length(amem_indx_min))"
end





