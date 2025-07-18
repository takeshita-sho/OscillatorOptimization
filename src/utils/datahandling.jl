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

