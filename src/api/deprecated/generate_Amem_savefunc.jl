function strip_t_from_symbol(sym)
    stripped_sym = Symbol(replace(string(sym), "(t)" => ""))
    return stripped_sym
end


function strip_t_from_string(str::String)
    stripped_str = replace(str, "(t)" => "")
    return stripped_str
end


function convert_symbolics_to_symbols(symbolics, sys)
    sys_name_map = get_name_to_var_dict(sys)
    sym_vec = [sys_name_map[s] for s in symbolics]
    return sym_vec
end


function replace_gamma_terms(str::String)
    replace(string(str), r"Γ\[(\d+)\]" => s"Γ\g<1>")
end


# Function to map eliminated variables (lhs of conserved_equations) to their right hand side (after zero_valued variables have been substituted)
function make_eliminated_var_mapping(conserved_equations, defaults)
    # Filter out zero-default variables
    zero_defaults = Dict(k => v for (k, v) in defaults if v isa Float64 && v == 0.0)

    # List to hold processed laws
    eliminated_var_mapping = Dict()

    for (i, eq) in enumerate(conserved_equations)
        # Substitute zero-default variables with 0 and simplify the expression
        # substituted_expr = substitute(eq.rhs, zero_defaults)
        substituted_expr = eq.rhs
        simplified_expr_str = string(simplify(substituted_expr))

        simplified_expr_str = strip_t_from_string(simplified_expr_str) |> replace_gamma_terms

        # Store the mapping
        stripped_t_lhs = string(strip_t_from_symbol(eq.lhs))
        eliminated_var_mapping[stripped_t_lhs] = simplified_expr_str
    end

    return eliminated_var_mapping
end


# Function to generate a dynamic function based on a system's symbols
function generate_Amem_expression_string(sys)
    # Retrieve system symbols and convert them to strings for filtering
    sys_syms = string.([get_sys_species_symbols(sys); get_observable_symbols(sys)])

    # Filter symbols for the numerator: include symbols containing "X" or "LpA"
    numerator_syms = filter(sym -> occursin("X", sym) || occursin("LpA", sym), sys_syms)

    # Filter symbols for the denominator: include symbols containing "A"
    denominator_syms = filter(sym -> occursin("A", sym) || occursin("X", sym), sys_syms)

    # Create string expressions for the sum of the numerator and denominator
    numerator_expr_str = join(numerator_syms, " + ")
    denominator_expr_str = join(denominator_syms, " + ")

    # Create the final expression string for Amem
    Amem_expr_str = "($numerator_expr_str) / ($denominator_expr_str)"

    return Amem_expr_str
end





# Function to replace variables in an expression string based on a mapping dictionary using regex for exact matches
function replace_eliminated_variables_in_expression(mapping::Dict{Any, Any}, expr_str::String)
    # Sort the keys by length in descending order to avoid partial replacement issues
    sorted_keys = sort(collect(keys(mapping)), by=length, rev=true)

    # Replace each variable in the expression string with its corresponding mapped value
    for key in sorted_keys
        # Create a regex pattern with word boundaries around the key
        pattern = Regex("\\b" * escape_string(key) * "\\b")
        println(pattern)

        # Check if the key occurs in the string and replace it
        # if occursin(pattern, expr_str)
        expr_str = replace(expr_str, pattern => "($(mapping[key]))")  # Use parentheses to ensure correct operation precedence
        # end
    end

    return expr_str
end


function add_X_multipliers(expression::String)
    # Regular expression to find terms containing 'X'
    pattern = r"(\b[X]{2,}\w*\b)"

    # Find all matches
    matches = collect(eachmatch(pattern, expression))

    # Create a dictionary to hold original terms and their multiplied versions
    replacements = Dict{String, String}()

    # Populate the dictionary with replacements
    for m in matches
        original = m.match
        x_count = count(c -> c == 'X', original)
        # Create the replacement string with the multiplier
        replacements[original] = "($(float(x_count)) * $original)"
    end

    # Perform the replacements
    multipled_expression = expression
    for (original, replacement) in replacements
        multipled_expression = replace(multipled_expression, Regex("\\b" * escape_string(original) * "\\b") => replacement)
    end

    return multipled_expression
end


function generate_subbed_Amem_expression_string(sys)
    Amem_expression_string = generate_Amem_expression_string(sys)

    conserved_eqs = get_observed(sys)
    sys_defaults = sys.defaults
    eliminated_var_mapping = make_eliminated_var_mapping(conserved_eqs, sys_defaults)

    Amem_expression_string = replace_eliminated_variables_in_expression(eliminated_var_mapping, Amem_expression_string)
    println(Amem_expression_string)

    X_multiplied_expression_string = add_X_multipliers(Amem_expression_string)
    println(X_multiplied_expression_string)

    # Convert to Symbolics expression in order to simplify 
    parsed_expr = Meta.parse(X_multiplied_expression_string)

    symbolic_expr = parse_expr_to_symbolic(parsed_expr, @__MODULE__)
    
    # Convert back to string
    Amem_expression_string = string(symbolic_expr)
end


# Generate the expression for (XXXBCD + XXBCD + XBCD + BCD) divided by all monomer species (B, C, D, and their intermediate complexes)
function generate_fulltrimer_expression_string(sys)
    # Retrieve system symbols and convert them to strings for filtering
    sys_syms = string.([get_sys_species_symbols(sys); get_observable_symbols(sys)])

    numerator_expr_str = "(XXXBCD + XXBCD + XBCD + BCD)" 

    # Filter symbols for the denominator: include symbols containing "B", "C", "D"
    denominator_syms = filter(sym -> occursin("B", sym) || occursin("C", sym) || occursin("D", sym), sys_syms)

    # Filter symbols for B, C, D complexes separately
    b_syms = filter(sym -> occursin("B", sym), sys_syms)
    c_syms = filter(sym -> occursin("C", sym), sys_syms)
    d_syms = filter(sym -> occursin("D", sym), sys_syms)

    # Create string expressions for the sum of each group
    b_expr_str = join(b_syms, " + ")
    c_expr_str = join(c_syms, " + ")
    d_expr_str = join(d_syms, " + ")

    # Create string expressions for the sum of the numerator and denominator
    # denominator_expr_str = join([b_expr_str, c_expr_str, d_expr_str], " + ")
    denominator_expr_str = "min($(b_expr_str), $(c_expr_str), $(d_expr_str))"

    # Create the final expression string for trimer assembly yield
    fulltrimer_expr_str = numerator_expr_str * " / " * "(" * denominator_expr_str * ")"

    return fulltrimer_expr_str
end

# # Generate the expression for (XXXBCD + XXBCD + XBCD + BCD) divided by the minimum of the sums of B, C, D species
# function generate_fulltrimer_expression_string(sys)
#     # Retrieve system symbols and convert them to strings for filtering
#     sys_syms = string.([get_sys_species_symbols(sys); get_observable_symbols(sys)])

#     numerator_expr_str = "(XXXBCD + XXBCD + XBCD + BCD)"

#     # Filter symbols for B, C, D complexes separately
#     b_syms = filter(sym -> occursin("B", sym), sys_syms)
#     c_syms = filter(sym -> occursin("C", sym), sys_syms)
#     d_syms = filter(sym -> occursin("D", sym), sys_syms)

#     # Create string expressions for the sum of each group
#     b_expr_str = join(b_syms, " + ")
#     c_expr_str = join(c_syms, " + ")
#     d_expr_str = join(d_syms, " + ")

#     # Use `min` to find the limiting group and form the denominator expression
#     denominator_expr_str = "min($(b_expr_str), $(c_expr_str), $(d_expr_str))"

#     # Create the final expression string for the yield
#     fulltrimer_expr_str = numerator_expr_str * " / " * "(" * denominator_expr_str * ")"

#     return fulltrimer_expr_str
# end



function generate_subbed_fulltrimer_expression_string(sys)
    fulltrimer_expression_string = generate_fulltrimer_expression_string(sys)
    println("Before substitution: ", fulltrimer_expression_string)

    conserved_eqs = get_observed(sys)
    sys_defaults = sys.defaults
    eliminated_var_mapping = make_eliminated_var_mapping(conserved_eqs, sys_defaults)

    fulltrimer_expression_string = replace_eliminated_variables_in_expression(eliminated_var_mapping, fulltrimer_expression_string)
    println(fulltrimer_expression_string)

    # Convert to Symbolics expression in order to simplify 
    # parsed_expr = Meta.parse(fulltrimer_expression_string)

    # symbolic_expr = parse_expr_to_symbolic(parsed_expr, @__MODULE__)
    
    # # Convert back to string
    # fulltrimer_expression_string = string(symbolic_expr)

    # println("After substitution: ",fulltrimer_expression_string)
    return fulltrimer_expression_string
end

# # Function to generate a dynamic function based on a list of system symbols
# function generate_Amem_callback(sys)
#     Amem_expression_string = generate_subbed_Amem_expression_string(sys)

#     # Convert to Symbolics expression in order to simplify 
#     Amem_subbed_expr = Meta.parse(Amem_expression_string)

#     Amem_symbolic_expr = parse_expr_to_symbolic(Amem_subbed_expr, @__MODULE__)

#     # Convert back to string
#     Amem_expression_string = string(Amem_symbolic_expr)

#     species_strings = string.(get_sys_species_symbols(sys))
#     parameter_strings = setdiff(string.(get_sys_symbols(sys)), species_strings)


#     # Sort the system strings by length in descending order to avoid partial replacement issues
#     sorted_species_strings = sort(species_strings, by=length, rev=true)
#     sorted_parameter_strings = sort(parameter_strings, by=length, rev=true)

#     # Create a dictionary to map each variable string to its index
#     # str_to_index = Dict(str => idx for (idx, str) in enumerate(sorted_sys_strings))
#     species_str_to_index = Dict(str => idx for (idx, str) in enumerate(species_strings))
#     parameter_str_to_index = Dict(str => idx for (idx, str) in enumerate(parameter_strings))


#     # Replace species symbols in the expression string with their indexed form
#     for str in sorted_species_strings
#         if occursin(str, Amem_expression_string)
#             idx = species_str_to_index[str]
#             Amem_expression_string = replace(Amem_expression_string, str => "u[$idx]")
#         end
#     end

#     for str in sorted_parameter_strings
#         if occursin(str, Amem_expression_string)
#             idx = parameter_str_to_index[str]
#             Amem_expression_string = replace(Amem_expression_string, str => "integrator.p[$idx]")
#         end
#     end

#     # Parse the modified expression string into a symbolic expression
#     Amem_expr = Meta.parse(Amem_expression_string)

#     # Construct the function expression
#     function_expr = quote
#         function save_func(u, t, integrator)
#             return $Amem_expr
#         end
#     end

#     return function_expr
# end



# Function to generate a dynamic function based on a list of system symbols
function generate_callback(expression_string, sys)

    println("Expression string: ", expression_string)


    species_strings = string.(get_sys_species_symbols(sys))
    println("Species strings: ", species_strings)
    parameter_strings = setdiff(string.(get_sys_symbols(sys)), species_strings)
    println("Parameter strings: ", parameter_strings)


    # Sort the system strings by length in descending order to avoid partial replacement issues
    sorted_species_strings = sort(species_strings, by=length, rev=true)
    sorted_parameter_strings = sort(parameter_strings, by=length, rev=true)
    
    # Replace all instances of "Γi" with "Γ[i]"
    sorted_parameter_strings = map(sorted_parameter_strings) do str
        if occursin("Γ", str)
            return replace(str, r"Γ(\d+)" => s"Γ[\1]")
        else
            return str
        end
    end
    
    println("Sorted parameter strings: ", sorted_parameter_strings)
    # Create a dictionary to map each variable string to its index
    # str_to_index = Dict(str => idx for (idx, str) in enumerate(sorted_sys_strings))
    species_str_to_index = Dict(str => idx for (idx, str) in enumerate(species_strings))
    # parameter_str_to_index = Dict(str => idx for (idx, str) in enumerate(parameter_strings))


    # Replace species symbols in the expression string with their indexed form
    for str in sorted_species_strings
        if occursin(str, expression_string)
            idx = species_str_to_index[str]
            expression_string = replace(expression_string, str => "u[$idx]")
        end
    end

    for str in sorted_parameter_strings
        if occursin(str, expression_string)
            # idx = parameter_str_to_index[str]
            expression_string = replace(expression_string, str => "integrator.ps[integrator.f.sys.$(str)]")
        end
    end

    # Parse the modified expression string into an Expression
    expr = Meta.parse(expression_string)

    # Construct the function expression
    function_expr = quote
        function save_func(u, t, integrator)
            return $expr
        end
    end

    return function_expr
end




