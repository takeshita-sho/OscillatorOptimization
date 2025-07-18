# Function to parse a block of lines and return a dictionary of molecule counts
function parse_time_block(block::Vector{String})
    # Define the dictionary mapping sorted tuples of (name, count) pairs to symbols
    molecule_to_symbol = Dict(
        # Single components
        (("lipid", 1),) => :Lp,
        (("kinase", 1),) => :K,
        (("phosphatase", 1),) => :P,
        (("adaptor", 1),) => :A,
        # Two components (tuples are sorted alphabetically by name)
        (("adaptor", 1), ("lipid", 1)) => :LpA,
        (("kinase", 1), ("lipid", 1)) => :LK,
        (("lipid", 1), ("phosphatase", 1)) => :LpP,
        (("adaptor", 1), ("kinase", 1)) => :AK,
        (("adaptor", 1), ("phosphatase", 1)) => :AP,
        # Three+ components (tuples are sorted alphabetically by name)
        (("adaptor", 1), ("kinase", 1), ("lipid", 1)) => :LpAK,
        (("adaptor", 1), ("lipid", 1), ("phosphatase", 1)) => :LpAP,
        (("adaptor", 1), ("kinase", 1), ("lipid", 2)) => :LpAKL, # Assuming "lipid: 2." means two lipids
        (("adaptor", 1), ("lipid", 2), ("phosphatase", 1)) => :LpAPLp # Assuming "lipid: 2." means two lipids
        # Add any other complex combinations here following the sorted tuple pattern
    )

    # Initialize a dictionary to store molecule counts for this time block
    molecule_counts = Dict{Symbol, Int}()

    # Regex to find "name: count." pairs, allowing for flexible spacing
    # Captures: 1=name (word), 2=count (digits)
    component_regex = r"(\w+):\s*(\d+)\."

    # Iterate over each line in the block, skipping the first line (time)
    for line in block[2:end]
        parts = split(line, '\t')
        if length(parts) != 2
            @warn "Skipping malformed line (expected 2 tab-separated parts): '$line'"
            continue
        end

        count_str, raw_complex_str = parts[1], strip(parts[2])

        count = tryparse(Int, count_str)
        if isnothing(count)
            @warn "Skipping line with non-integer count: '$line'"
            continue
        end

        if isempty(raw_complex_str)
            # Handle case where complex string might be empty for some reason
            continue
        end

        # --- Canonical Key Generation (Tuple of Pairs) ---
        components = Tuple{String, Int}[]
        last_match_end = 0
        parse_error = false
        for match in eachmatch(component_regex, raw_complex_str)
            # Check if matches are contiguous (accounts for spaces between components)
            if match.offset > last_match_end + 1 && last_match_end != 0
                 # Check if the gap is just whitespace
                 gap = SubString(raw_complex_str, last_match_end + 1, match.offset - 1)
                 if !all(isspace, gap)
                    @warn "Unexpected characters between components in '$raw_complex_str'. Skipping line: '$line'"
                    parse_error = true
                    break
                 end
            elseif match.offset == 0 && last_match_end !=0 # Should not happen with correct regex use but safety check
                 @warn "Regex match overlap detected in '$raw_complex_str'. Skipping line: '$line'"
                 parse_error = true
                 break
            end


            component_name = match.captures[1]
            component_count_str = match.captures[2]
            component_count = parse(Int, component_count_str) # Should not fail due to regex digit match

            push!(components, (component_name, component_count))
            last_match_end = match.offset + length(match.match) -1
        end

        # Check if the entire string was consumed by matches + whitespace
        if !parse_error && last_match_end < length(raw_complex_str)
            remaining = SubString(raw_complex_str, last_match_end + 1)
             if !all(isspace, remaining) # Allow trailing whitespace
                 @warn "Trailing unexpected characters found in '$raw_complex_str' after parsing components. Skipping line: '$line'"
                parse_error = true
            end
        end

        if parse_error || isempty(components) # Skip if errors or no components found
            # Check if raw_complex_str was not empty but parsing failed
            if !parse_error && !isempty(raw_complex_str)
                 @warn "Could not parse any components from non-empty string: '$raw_complex_str' in line '$line'"
            end
            continue
        end


        # Sort the components alphabetically by name
        sort!(components, by = first)

        # Convert the sorted list of pairs into an immutable Tuple
        key_tuple = Tuple(components)
        # --- End Canonical Key Generation ---

        # Find the corresponding symbol using the canonical tuple key
        symbol = get(molecule_to_symbol, key_tuple, :Unknown)

        if symbol != :Unknown
            molecule_counts[symbol] = get(molecule_counts, symbol, 0) + count
        else
            # Log unknown complexes for debugging if needed
            @warn "Unknown complex encountered. Key: $key_tuple (Raw: '$raw_complex_str')"
        end
    end

    return molecule_counts
end

# Example Usage (assuming `block` is a Vector{String} like in the file):
# test_block = [
#     "Time (s): 0.065",
#     "429\tlipid: 1.",
#     "466\tadaptor: 1. lipid: 1.", # Order 1
#     "133\tadaptor: 1.",
#     "236\tphosphatase: 1.",
#     "299\tkinase: 1.",
#     "2\tlipid: 1. adaptor: 1. kinase: 1.", # Order 2
#     "4\tphosphatase: 1. lipid: 1.",
#     "1\tadaptor: 1. phosphatase: 1. lipid: 2."
# ]
# counts = parse_time_block(test_block)
# println(counts)

# Function to extract timepoints from the file
function extract_timepoints(filepath::String)
    timepoints = Float64[]
    open(filepath, "r") do file
        for line in eachline(file)
            if startswith(line, "Time (s):")
                time = parse(Float64, split(line, ':')[2])
                push!(timepoints, time)
            end
        end
    end
    return timepoints
end

# Function to initialize the DataFrame
function initialize_dataframe(filepath::String)
    # Extract timepoints from the file
    timepoints = extract_timepoints(filepath)
    
    # Define the columns of the DataFrame
    columns = [:Time, :Lp, :K, :P, :A, :LpA, :LK, :LpP, :LpAK, :LpAP, :LpAKL, :LpAPLp, :AK, :AP]
    
    # Create a DataFrame with the correct size, initializing non-Time columns to zero
    df = DataFrame()
    df.Time = timepoints
    for col in columns[2:end]  # Skip the Time column
        df[!, col] .= 0
    end
    
    return df
end



# Function to read the file and return an array of time blocks
function extract_time_blocks(filepath::String)
    time_blocks = Vector{Vector{String}}()  # Array of arrays to hold blocks of lines
    current_block = Vector{String}()        # Temporary storage for the current block of lines

    open(filepath, "r") do file
        for line in eachline(file)
            if startswith(line, "Time (s):")  # Check if the line indicates a new time block
                if !isempty(current_block)      # If the current block is not empty, save it
                    push!(time_blocks, current_block)
                    current_block = Vector{String}()  # Start a new block
                end
            # molecule_string = split(line, '\t')[2]
        end
        push!(current_block, line)  # Add line to the current block
        end
        if !isempty(current_block)  # Add the last block if it's not empty
            push!(time_blocks, current_block)
        end
    end

    return time_blocks
end



# Main function to process the file and update the DataFrame
function process_and_update_dataframe(df::DataFrame, time_blocks::Vector{Vector{String}})
    # Ensure that the number of time blocks matches the number of rows in the DataFrame
    if size(df, 1) != length(time_blocks)
        throw(ArgumentError("The number of DataFrame rows and time blocks must match"))
    end

    # Loop through each time block and corresponding DataFrame row
    for (index, block) in enumerate(time_blocks)
        # Parse the current time block to get molecule counts
        molecule_counts = parse_time_block(block)

        # Update the DataFrame row with the molecule counts
        for (molecule, count) in molecule_counts
            df[index, molecule] = count
        end
    end
end


function compute_Amem(df::AbstractDataFrame)
    # Use Regex to select columns that contain "LpA" in their name
    LpA_columns = select(df, r"LpA")
    sum_LpA = sum(eachcol(LpA_columns))

    # Use Regex to select columns that contain "A" in their name
    A_columns = select(df, r"A")
    sum_A = sum(eachcol(A_columns))

    # Compute Amem
    Amem = sum_LpA ./ sum_A

    return Amem
end



# ============================================================================
#  Fast histogram parser: column-wise construction (no per-cell mutation)
# ============================================================================

"""
    make_hist_dataframe(path::String) -> DataFrame

Parses the histogram_complexes_time.dat file and returns a DataFrame with the timepoints and the counts of each complex.
"""
function make_hist_dataframe(path::String)
    @assert occursin("histogram", path)

    _SYMBOLS = [:Lp, :K, :P, :A, :LpA, :LK, :LpP, :LpAK, :LpAP, :LpAKL, :LpAPLp, :AK, :AP]

    # single streaming pass: build vectors incrementally
    times = Float64[]
    col_data = Dict{Symbol, Vector{Int}}()
    for s in _SYMBOLS
        col_data[s] = Int[]
    end
    current_counts = Dict{Symbol,Int}()
    last_time = 0.0
    open(path) do io
        for ln in eachline(io)
            if startswith(ln, "Time (s):")
                # dump previous counts (first block skip)
                if !isempty(current_counts)
                    push!(times, last_time)
                    for s in _SYMBOLS
                        push!(col_data[s], get(current_counts, s, 0))
                    end
                    empty!(current_counts)
                end
                last_time = parse(Float64, ln[11:end])
            elseif !isempty(ln)
                sp = split(ln, '\t'; limit = 2)
                length(sp) == 2 || continue
                count = parse(Int, sp[1])
                raw   = sp[2]

                # Fast order-independent detection via bitmask
                has_lip1 = occursin("lipid: 1.", raw)
                has_lip2 = occursin("lipid: 2.", raw)
                has_ada  = occursin("adaptor: 1.", raw)
                has_kin  = occursin("kinase: 1.",  raw)
                has_pho  = occursin("phosphatase: 1.", raw)

                sym = nothing
                if has_lip1 && !(has_ada || has_kin || has_pho)
                    sym = :Lp
                elseif has_kin && !(has_ada || has_lip1 || has_pho)
                    sym = :K
                elseif has_pho && !(has_ada || has_lip1 || has_kin)
                    sym = :P
                elseif has_ada && !(has_lip1 || has_kin || has_pho)
                    sym = :A
                elseif has_ada && has_lip1 && !(has_kin || has_pho)
                    sym = :LpA
                elseif has_kin && has_lip1 && !(has_ada || has_pho)
                    sym = :LK
                elseif has_pho && has_lip1 && !(has_ada || has_kin)
                    sym = :LpP
                elseif has_ada && has_kin && has_lip1 && !(has_pho)
                    sym = :LpAK
                elseif has_ada && has_pho && has_lip1 && !(has_kin)
                    sym = :LpAP
                elseif has_ada && has_kin && has_lip2 && !(has_pho)
                    sym = :LpAKL
                elseif has_ada && has_pho && has_lip2 && !(has_kin)
                    sym = :LpAPLp
                elseif has_ada && has_kin && !(has_lip1 || has_lip2 || has_pho)
                    sym = :AK
                elseif has_ada && has_pho && !(has_lip1 || has_lip2 || has_kin)
                    sym = :AP
                end

                if sym !== nothing
                    current_counts[sym] = get(current_counts, sym, 0) + count
                end
            end
        end
    end
    # dump last block
    if !isempty(current_counts)
        push!(times, last_time)
        for s in _SYMBOLS
            push!(col_data[s], get(current_counts, s, 0))
        end
    end

    cols = (; Time = times, (s => col_data[s] for s in _SYMBOLS)...)
    df = DataFrame(cols; copycols=false)
    df.Amem = compute_Amem(df)
    return df
end
