# using Makie
using LaTeXStrings
"""
    parse_parameter(param::String) -> NamedTuple

Parse a parameter name (rate constant or species) and return its formatted properties.
Example: "kfᴸᴷ" -> (formatted_name="kf_LK", color=:green, unit="μM⁻¹ s⁻¹")
        "L" -> (formatted_name="L", color=:orange, unit="μM")

# Arguments
- `param`: Raw parameter name (e.g., "kfᴸᴷ", "krᴬᴾ", "kcatᴸᴾ", "L", "K", "P", "A")
"""
function parse_parameter(param::String)
    # Dictionary definitions
    species_colors = Dict("L" => :orange, "K" => :teal, "P" => :purple, "A" => :red)
    pattern_colors = Dict("ᴸᴬ" => :red, "ᴸᴷ" => :teal, "ᴸᴾ" => :purple, "ᴬᴷ" => :pink, "ᴬᴾ" => :green)
    units = Dict(
        "kf" => "μM⁻¹ s⁻¹", 
        "kr" => "s⁻¹", 
        "kcat" => "s⁻¹", 
        "KD" => "μM",
        "KM" => "μM",
        "kE" => "μM⁻¹ s⁻¹",
        "Amem" => "dimensionless"
    )
    
    # Handle special case for Amem
    if param == "Amem"
        return (
            formatted_name = L"$\mathit{A}_{\mathrm{mem}}$",
            color = species_colors["A"],  # Use the color for A species
            unit = units["Amem"]
        )
    end
    
    # Handle species case
    if param in ["L", "K", "P", "A"]
        return (
            formatted_name = L"$\mathrm{%$param}$",
            color = species_colors[param],
            unit = "μM"
        )
    end
    
    # Extract pattern and handle rate constants
    pattern = match(r"(ᴸᴬ|ᴸᴾ|ᴸᴷ|ᴬᴷ|ᴬᴾ)", param)
    isnothing(pattern) && error("Unrecognized parameter pattern in: $param")
    pattern = String(pattern.match)
    
    # Convert pattern to standard notation (not unicode superscript)
    # Create a mapping dictionary for small capitals to regular letters
    small_caps_map = Dict('ᴸ' => 'L', 'ᴬ' => 'A', 'ᴷ' => 'K', 'ᴾ' => 'P')
    label = join([get(small_caps_map, c, c) for c in pattern])
    
    # Debug print
    # println("Parameter: $param, Pattern: $pattern, Label: $label")
    
    # Get parameter type and unit
    param_type = if startswith(param, "k")
        String(match(r"^(kf|kr|kcat)", param).match)
    elseif startswith(param, "Kd")
        "KD"
    elseif startswith(param, "Km")
        "KM"
    elseif startswith(param, "Ke")
        "kE"
    else
        error("Unrecognized parameter type: $param")
    end
    
    # Format name with LaTeX - explicitly set variables to italic and use solidus for fractions
    formatted_name = if param_type == "kE"
        # Catalytic efficiency: k_cat/K_M with solidus
        # Use a small negative space (\!) on one side only for better spacing
        LaTeXString("\\mathit{k}_{\\mathrm{cat}}^{$label}/\\mathit{K}_{\\mathrm{M}}^{$label}")
    elseif param_type == "KD"
        LaTeXString("\\mathit{K}_{\\mathrm{D}}^{$label}")
    elseif param_type == "KM"
        LaTeXString("\\mathit{K}_{\\mathrm{M}}^{$label}")
    elseif param_type == "kf"
        LaTeXString("\\mathit{k}_{\\mathrm{f}}^{$label}")
    elseif param_type == "kr"
        LaTeXString("\\mathit{k}_{\\mathrm{r}}^{$label}")
    elseif param_type == "kcat"
        LaTeXString("\\mathit{k}_{\\mathrm{cat}}^{$label}")
    else
        error("Unrecognized parameter type: $param_type")
    end
    
    return (
        formatted_name = formatted_name,
        color = pattern_colors[pattern],
        unit = units[param_type]
    )
end


# function test_parameter_parsing()
#     # Test cases in a logical grouping
#     test_params = [
#         # Group 1: LK pathway
#         "kfᴸᴷ", "krᴸᴷ", "kcatᴸᴷ", "Kmᴸᴷ", "Keᴸᴷ",
#         # Group 2: LP pathway
#         "kfᴸᴾ", "krᴸᴾ", "kcatᴸᴾ", "Kmᴸᴾ", "Keᴸᴾ",
#         # Group 3: LA and other patterns
#         "kfᴸᴬ", "krᴸᴬ", "Kdᴸᴬ",
#         "kfᴬᴷ", "krᴬᴷ", "Kdᴬᴷ",
#         "kfᴬᴾ", "krᴬᴾ", "Kdᴬᴾ",
#         # Group 4: Species
#         "L", "K", "P", "A", "Amem"
#     ]
    
#     # Create figure
#     fig = Figure(size=(800, 400))
#     ax = Axis(fig[1,1];
#         title = "Parameter Formatting Test",
#         xticklabelsvisible = false,
#         yticklabelsvisible = false,
#         xgridvisible = false,
#         ygridvisible = false
#     )
    
#     # Calculate grid positions
#     n_rows = 6
#     n_cols = 4
#     x = repeat(1:n_cols, n_rows)
#     y = repeat(n_rows:-1:1, inner=n_cols)
    
#     # Get formatted names
#     formatted_names = [parse_parameter(param).formatted_name for param in test_params]
    
#     # Plot text in grid
#     text!(ax, x[1:length(test_params)], y[1:length(test_params)];
#         text = formatted_names,
#         align = (:center, :center),
#         fontsize = 25
#     )
    
#     # Set axis limits with padding
#     xlims!(ax, 0.5, n_cols + 0.5)
#     ylims!(ax, 0.5, n_rows + 0.5)
    
#     hidespines!(ax)
    
#     return fig
# end

# using CairoMakie
# test_parameter_parsing()