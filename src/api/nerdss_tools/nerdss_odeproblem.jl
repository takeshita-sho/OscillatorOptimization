# Extract WaterBox dimensions
function extract_waterbox(file_content::String)
    
    # Find the line containing "WaterBox"
    waterbox_line = match(r"WaterBox\s*=\s*\[(\d+),\s*(\d+),\s*(\d+)\]", file_content)
    
    if waterbox_line === nothing
        error("WaterBox line not found in the file.")
    end
    
    waterbox_dims = [parse(Int, dim) for dim in waterbox_line.captures[1:3]]
    
    return Float64.(waterbox_dims)
end


function extract_rate_constants(file_content::String)
    
    # Regex to find lines with rate constants marked by comments starting with "#k"
    pattern = r"(\w+)\s*=\s*([\d\.]+)\s*#\s*(k\w+)"
    matches = eachmatch(pattern, file_content)

    new_k_table = Dict(
        "ka1" => :kfᴸᴷ,
        "kb1" => :krᴸᴷ,
        "kcat1" => :kcatᴸᴷ,
        "ka2" => :kfᴸᴬ,
        "kb2" => :krᴸᴬ,
        "ka3" => :kfᴬᴷ,
        "kb3" => :krᴬᴷ,
        "ka4" => :kfᴬᴾ,
        "kb4" => :krᴬᴾ,
        "ka7" => :kfᴸᴾ,
        "kb7" => :krᴸᴾ,
        "kcat7" => :kcatᴸᴾ,
    )
    
    # Extract rate constants and their values, converting to appropriate types
    rate_constants = [get(new_k_table, string(m.captures[3]), Symbol(m.captures[3])) => parse(Float64, m.captures[2]) for m in matches]
    
    return rate_constants
end


function rate_converter(conc::Float64, volume_nm3::Float64)
    # Convert concentration from per microsecond to per second
    conc_per_sec = conc * 1e6
    
    # # Convert volume from cubic micrometers (um^3) to cubic nanometers (nm^3)
    # volume_nm3 = volume_um3 * 1e9
    
    # Calculate and return the rate in copies per second
    return conc_per_sec / volume_nm3
end


function convert_forward_rate_constants(rate_constants::Vector{Pair{Symbol, Float64}}, volume_nm3::Float64)
    # Initialize an empty vector to store the converted rate constants
    converted_rate_constants = Vector{Pair{Symbol, Float64}}()
    # converted_rate_constants = similar(rate_constants)

    # Iterate over each rate constant
    for (rate_constant_symbol, rate_constant_value) in rate_constants
        # Check if the symbol starts with "kf"
        if startswith(string(rate_constant_symbol), "kf")
            # Convert the rate constant value
            converted_value = rate_converter(rate_constant_value, volume_nm3)
            # Append the converted rate constant to the vector
            push!(converted_rate_constants, rate_constant_symbol => converted_value)
        else
            # Append the unmodified rate constant to the vector
            push!(converted_rate_constants, rate_constant_symbol => rate_constant_value)
        end
    end

    return converted_rate_constants
end




function extract_initial_copynumbers(file_content::String)
    # Extract the block of text between "start molecules" and "end molecules"
    molecules_block = match(Regex("start molecules\\s*(.*?)\\s*end molecules", "ms"), file_content)
    
    if molecules_block === nothing
        error("Molecules block not found in the file.")
    end
    
    # Extract individual molecule lines
    lines = split(strip(molecules_block.captures[1]), '\n')
    
    # Initialize an empty vector for the results
    molecule_counts = Vector{Pair{Symbol, Int}}()
    
    # Process each line to extract molecule counts and map them to the specified symbols
    for line in lines
        if occursin("pip2", line)
            counts = match(r"pip2\s*:\s*(\d+)\s*\(head~U\)\s*(\d+)\s*\(head~P\)", line)
            push!(molecule_counts, :L => parse(Int, counts.captures[1]))
            push!(molecule_counts, :Lp => parse(Int, counts.captures[2]))
        elseif occursin("ap2", line)
            count = match(r"ap2\s*:\s*(\d+)", line)
            push!(molecule_counts, :A => parse(Int, count.captures[1]))
        elseif occursin("kin", line)
            count = match(r"kin\s*:\s*(\d+)", line)
            push!(molecule_counts, :K => parse(Int, count.captures[1]))
        elseif occursin("syn", line)
            count = match(r"syn\s*:\s*(\d+)", line)
            push!(molecule_counts, :P => parse(Int, count.captures[1]))
        end
    end
    
    return molecule_counts
end


# Function to extract simulation parameters and compute the total simulation time
function extract_simulation_time(file_content::String)
    # Extract nItr and timeStep using regular expressions
    nItr_match = match(r"nItr\s*=\s*(\d+)", file_content)
    timeStep_match = match(r"timeStep\s*=\s*(\d+)", file_content)
    
    if nItr_match === nothing || timeStep_match === nothing
        error("Simulation parameters not found in the file.")
    end
    
    # Parse the extracted values
    nItr = parse(Int, nItr_match.captures[1])
    timeStep = parse(Int, timeStep_match.captures[1])
    
    # Calculate the total time in seconds
    total_time_seconds = nItr * timeStep / 1e6
    
    return (0.0, total_time_seconds)
end



#< Main function to setup the ODE problem
function make_nerdss_ode_problem(directory_path::String)
    # Construct the full path to the parameters file
    parms_path = joinpath(directory_path, "parms.inp")

    # Read the entire content of the file once
    file_content = Base.read(parms_path, String)

    # Extract WaterBox dimensions
    waterbox_dims = extract_waterbox(file_content)
    volume = prod(waterbox_dims) # Calculate volume in nm^3

    # Extract rate constants
    rate_constants = extract_rate_constants(file_content)
    push!(rate_constants, :DF => waterbox_dims[end]) # Append DF to rate constants

    # Convert rate constants based on the volume
    converted_rate_constants = convert_forward_rate_constants(rate_constants, volume)

    # Extract initial molecule counts
    initial_copynumbers = extract_initial_copynumbers(file_content)

    # Extract simulation time span
    tspan = extract_simulation_time(file_content)

    # Setup the ODE problem
    nerdss_odeprob = OscTools.make_odeprob(;rn = fullrn, tspan = tspan)

    nerdss_odeprob = remake(nerdss_odeprob, u0 = initial_copynumbers, p = converted_rate_constants)
    # replace!(Tunable(), nerdss_odeprob.u0, initial_copynumbers)
    # replace!(Tunable(), parameter_values(nerdss_odeprob), converted_rate_constants)

    return nerdss_odeprob
end

