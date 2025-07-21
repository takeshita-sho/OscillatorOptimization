# Unit test for the fitness function calculation
# Tests that the fitness function produces non-zero results for a known oscillatory solution

using Test, Random
using OscillatorOptimization
using OrdinaryDiffEq: Rodas5, solve
using FFTW

@testset "Fitness function test" begin
    # Use the same solution from solveODE_test that we know is oscillatory
    sol = solve(make_odeprob(), Rodas5(), saveat = 0.1, abstol = 1e-7, reltol = 1e-7, save_on = true)
    
    # Extract the observables we need - both Amem_old and Amem
    Amem_old = sol[:Amem_old]  # For FFT analysis
    Amem = sol[:Amem]          # For time-domain analysis  
    time_vec = sol.t
    
    @test argmax(Amem) == 25  # Verify we have the same solution as solveODE_test
    
    # Create the saved_array structure that fitness function expects (2×N)
    saved_array = [Amem_old'; Amem']  # 2×length(time_vec) array
    
    @test size(saved_array) == (2, length(time_vec))
    
    # Create minimal components needed for fitness calculation
    # Create FFTW plan matching what saved_array[1, :] will produce
    sample_fft_array = @view saved_array[1, :]
    rfft_plan = plan_rfft(sample_fft_array)
    
    # Create minimal optsys-like object with just what we need
    optsys_minimal = (rfft_plan = rfft_plan, t = time_vec)
    
    # Test the fitness function
    phenotype = zeros(3)
    result = calculate_fitness!(phenotype, saved_array, optsys_minimal)
    
    @test result === phenotype  # Function should return the same array
    
    # Test for specific expected values from known oscillatory solution
    expected_fitness = 0.520
    expected_period = 49.2
    expected_amplitude = 0.630
    
    println("Fitness function test results:")
    println("  Expected: Fitness=$(expected_fitness), Period=$(expected_period), Amplitude=$(expected_amplitude)")
    println("  Actual:   Fitness=$(round(phenotype[1], digits=3)), Period=$(round(phenotype[2], digits=1)), Amplitude=$(round(phenotype[3], digits=3))")
    
    # Test for approximate equality with reasonable tolerance
    @test isapprox(phenotype[1], expected_fitness, rtol=0.01)  # 1% tolerance for fitness
    @test isapprox(phenotype[2], expected_period, rtol=0.01)   # 1% tolerance for period  
    @test isapprox(phenotype[3], expected_amplitude, rtol=0.01) # 1% tolerance for amplitude
    
    println("  ✓ All fitness function values match expected results within 1% tolerance")
end