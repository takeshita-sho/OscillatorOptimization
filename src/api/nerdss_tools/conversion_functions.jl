function convert_forward_rate_bounds(bounds_dict::Dict{Symbol, Tuple{Float64, Float64}}, volume_nm3::Float64)
    converted_bounds_dict = deepcopy(bounds_dict)

    # Iterate over each rate constant
    for (rate_constant_symbol, rate_constant_bounds) in bounds_dict
        # Check if the symbol starts with "kf"
        if startswith(string(rate_constant_symbol), "kf")
            # Convert the rate constant value
            converted_bounds = rate_converter.(rate_constant_bounds, volume_nm3)
            # Append the converted rate constant to the vector
            converted_bounds_dict[rate_constant_symbol] = converted_bounds
        else
            # Append the unmodified rate constant to the vector
            converted_bounds_dict[rate_constant_symbol] = rate_constant_bounds
        end
    end

    return converted_bounds_dict
end



function concentration_to_copy_number(concentration_μM, volume_nm3)
    # Avogadro's number (in molecules/mol)
    N_A = 6.022e23
    
    # Convert volume from nm^3 to liters
    # 1 nm^3 = 1e-24 liters
    volume_L = volume_nm3 * 1e-24
    
    # Convert concentration from μM to M
    # 1 μM = 1e-6 M
    concentration_M = concentration_μM * 1e-6
    
    # Calculate the copy number: C = N_A * V * c
    copy_number = N_A * volume_L * concentration_M
    
    # Return the integer value of the copy number
    return round(Int, copy_number) # Round and convert to integer for a discrete number of molecules
end

function copy_number_to_concentration(copy_number::Int, volume_nm3::Float64)
    # Avogadro's number (in molecules/mol)
    N_A = 6.022e23
    
    # Convert copy numbers to micromoles 
    micromoles = (copy_number / N_A) * 1e-6

    # Convert volume from nm^3 to liters
    volume_L = volume_nm3 * 1e-24

    # Calculate concentration in μM
    concentration_μM = micromoles / volume_L
    return concentration_μM
end






function convert_species_bounds(bounds_dict::Dict{Symbol, Tuple{Float64, Float64}}, volume_nm3::Float64)
    converted_bounds_dict = deepcopy(bounds_dict)

    for (species_symbol, species_bounds) in bounds_dict
        if !startswith(string(species_symbol), "k")
            converted_bounds_dict[species_symbol] = concentration_to_copy_number.(species_bounds, volume_nm3)
        else
            converted_bounds_dict[species_symbol] = species_bounds
        end
    end

    return converted_bounds_dict
end