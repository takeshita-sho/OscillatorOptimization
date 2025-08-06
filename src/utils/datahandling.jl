"""
    read_all_csvs(sweepDF_datapath)

Reads in all the csvs of all replicates in the given directory and concatenates them into a single dataframe.
"""
function read_all_csvs(sweepDF_datapath)
    # Get all replicate directories
    replicate_dirs = filter(x -> startswith(basename(x), "replicate_"), readdir(sweepDF_datapath, join=true))
    
    # Pre-allocate vector to store DataFrames
    dfs = Vector{DataFrame}()
    
    for dir in replicate_dirs
        # Extract replicate number from directory name
        replicate_num = parse(Int, match(r"replicate_(\d+)", basename(dir))[1])
        
        # Get all CSV files in this replicate directory
        csv_files = filter(x -> endswith(x, ".csv"), readdir(dir, join=true))
        
        # Read all CSVs in this replicate directory at once
        if !isempty(csv_files)
            df = CSV.read(csv_files, DataFrame)
            df.replicate .= replicate_num  # Broadcast assignment is more efficient
            push!(dfs, df)
        end
    end
    
    # Combine all DataFrames vertically using reduce
    return reduce(vcat, dfs)
end


"""
    solve_row(row::DataFrameRow, odeprob::ODEProblem; tspan = (0.0, 2000.0), abstol = 1e-10, reltol = 1e-10) -> ODESolution

Solves a single row of a dataframe using the given ODE problem. For interactive use; the regex makes it not performant.
"""
function solve_row(row::DataFrameRow, odeprob::ODEProblem; tspan = (0.0, 2000.0), abstol = 1e-10, reltol = 1e-10, kwargs...)
    # Get the parameters from the row 
    params = row[r"k|DF"] |> pairs |> collect

    # Get the initial conditions from the row 
    init_conds = row[r"^(L|K|P|A|B|C|D)$"] |> pairs |> collect
    
    # Remake the problem 
    newprob = remake(odeprob, p = params, u0 = init_conds)

    # Solve the problem 
    sol = solve(newprob, Rodas5P(autodiff = AutoForwardDiff(chunksize = length(odeprob.u0))), abstol = abstol, reltol = reltol, tspan = tspan, kwargs...)

    return sol
end

"""
    solve_row(data::DataFrame, odeprob::ODEProblem; tspan = (0.0, 2000.0), abstol = 1e-10, reltol = 1e-10, kwargs...) -> Function

Returns a performant function that can solve ODEs for any row index in the DataFrame using SymbolicIndexingInterface.
The returned function takes a row index and returns an ODESolution.

# Arguments
- `data`: DataFrame containing parameters and initial conditions
- `odeprob`: Base ODE problem to modify

# Returns
- `Function`: A function `f(row_index)` that returns an ODESolution for the specified row

# Example
```julia
# Create the solver function
solver_func = solve_row(df, odeprob)

# Solve for specific rows
sol1 = solver_func(1)  # Solve first row
sol42 = solver_func(42)  # Solve 42nd row
```
"""
function solve_row(data::AbstractDataFrame, odeprob::ODEProblem; tspan = (0.0, 2000.0), abstol = 1e-10, reltol = 1e-10, kwargs...)
    # Extract columns using existing regex (only once)
    parameter_cols = r"k|DF"
    species_cols = r"^(L|K|P|A|B|C|D)$"

    # Pre-extract data to matrices for fast indexing
    param_matrix = Matrix(data[!, parameter_cols])
    species_matrix = Matrix(data[!, species_cols])

    # Create fast setters
    set_params = setp(odeprob, parameter_cols)
    set_species = setu(odeprob, species_cols)
    
    # Return function that uses fast setters
    return function(row_index::Int; kwargs...)
        
        # Create thread-safe problem copy
        prob = remake(odeprob, p = copy(odeprob.p), u0 = copy(odeprob.u0))
        
        # Set the parameters and species using fast matrix indexing
        set_params(prob, view(param_matrix, row_index, :))
        set_species(prob, view(species_matrix, row_index, :))
        
        # Solve and return
        return solve(prob, Rodas5P(autodiff = AutoForwardDiff(chunksize = length(odeprob.u0))), abstol = abstol, reltol = reltol, tspan = tspan, kwargs...)
    end
end


