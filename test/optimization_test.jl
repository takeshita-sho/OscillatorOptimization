# Tests that standard optimization runs keep producing the same results

using Random, Test
using OscillatorOptimization
using DataFrames: nrow

@testset "Optimization test" begin
    fixed_params = Dict{Symbol, Float64}(:DF => 20.0)
    opt_sys = OptimizationReactionSystem(fullrn; fixed_params)
    results = run_optimization(10000, opt_sys;
        mutationRate = 0.95,     # Probability of applying mutation to an individual (0-1)
        crossoverRate = 0.75,    # Probability of applying crossover between parents (0-1)
        mutation_δ = 1.0,        # Mutation step size - larger values = bigger mutations
        pm = 0.25,              # Probability of mutating each gene
        η = 1,                  # Distribution index - lower values = more diverse mutations
        sbx_pm = 0.3,           # Probability of gene-wise mutation during crossover
        sbx_η = 1,              # Distribution index for crossover - lower = more diverse offspring
        num_tournament_groups = 20,  # Number of tournament groups - affects selection pressure
        n_points = Inf,         # Stop when this many unique oscillatory solutions found

        seed = 1234,           # Random seed for reproducibility, uses MersenneTwister
        show_trace = false,     # Show optimization progress
        parallelization = :thread  # Use multi-threading for evaluations and show progress bar in stdout
    )

    # Check that the number of unique oscillatory solutions is correct
    @test nrow(results.df) == 643

    @test results.df[end, :Fitness] ≈ 0.116003 rtol=1e-5
    @test results.df[end, :Period] ≈ 695.0 rtol=1e-5
end 