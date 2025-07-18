# Function to dynamically generate the conservation law computation function for DimensionalData structures
function generate_conservation_law_function(reaction_system::ReactionSystem)
    # Get the conservation laws from the reaction system
    conserved_quantities = conservationlaw_constants(reaction_system) # This returns an expression for each conserved quantity, which are elements of a Symbolics array and indexed like Γ[1], Γ[2], etc.

    # Because only some species ever have non-zero initial values, this just removes those zero-valued species from the conserved quantity expressions, and simplifies the expressions
    processed_conserved_quantities = preprocess_conservation_laws(conserved_quantities, reaction_system.defaults)
    
    # Start building the function expression
    function_expr = Expr(:function, Expr(:call, :compute_constants!, :individual))
    body_expr = Expr(:block)
    
    # Iterate over each conservation law and build the function body
    for law in processed_conserved_quantities
        var = law[1]
        expr = law[2]  # Right-hand side of the equation
        
        # Process the expression to handle operators and variables correctly
        processed_expr = process_expression(expr)
        
        # Use the Genes dimension with the At() selector for indexing into `individual`
        assignment_expr = :($(Expr(:ref, :individual, Expr(:kw, :Genes, :(At($(QuoteNode(var))))))) .= $processed_expr)
        
        push!(body_expr.args, assignment_expr)
    end
    
    # Add a return statement
    push!(body_expr.args, :nothing)
    
    # Complete the function expression
    push!(function_expr.args, body_expr)
    
    # Print the function expression for debugging
    println("Generated function expression:")
    println(function_expr)
    
    # Return the entire function as an expression
    return function_expr
end

# Helper function to process expressions, stripping time dependencies and formatting for DimensionalData
function process_expression(expr)
    expr_str = replace(string(expr), "(t)" => "")
    new_expr = Meta.parse(expr_str)
    return replace_variables_with_at(new_expr)
end


# Recursively replace variables in an expression with their Genes dimension indexed form, and replace operators with their dotted forms for broadcasting
function replace_variables_with_at(expr)
    if expr isa Symbol
        # Check if the symbol is an operator, if so, convert it to its broadcasting form
        if expr in (:+, :-, :*, :/, :^)
            # return Symbol('.' * string(expr))  # Convert to dotted operator
            return expr
        else
            return Expr(:ref, :individual, Expr(:kw, :Genes, :(At($(QuoteNode(expr))))))
        end
    elseif expr isa Expr
        # Apply the function recursively to all arguments of the expression
        expr.args = replace_variables_with_at.(expr.args)
        return expr
    else
        return expr
    end
end




# Function to preprocess conservation laws using Symbolics.jl
function preprocess_conservation_laws(conservation_laws, defaults)
    # Filter out zero-default variables
    zero_defaults = Dict(k => v for (k, v) in defaults if v == 0.0)

    # List to hold processed laws
    processed_laws = []

    for (i, law) in enumerate(conservation_laws)
        # Substitute zero-default variables with 0 and simplify the expression
        substituted_expr = substitute(law.rhs, zero_defaults)
        simplified_expr = simplify(substituted_expr)
        
        # Store the processed law with a new variable symbol
        push!(processed_laws, (Symbol("Γ", i), simplified_expr))
    end

    return processed_laws
end