#!/usr/bin/env julia

using OscillatorOptimization
using Random

println("Julia version: ", VERSION)

# Quick test - create a single individual and check fitness
fixed_params = Dict{Symbol, Float64}(:DF => 20.0)
opt_sys = OptimizationReactionSystem(fullrn; fixed_params)

rng = MersenneTwister(1234)
test_population = generate_population(fullrn, 1, fixed_params; rng=rng)
test_individual = first(test_population)  # SlicedDimArray acts like a 1D vector

println("Testing individual evaluation...")
sol, saved_array = solve_odes(test_individual, opt_sys)
println("ODE solution status: ", sol.retcode)
println("Saved array size: ", size(saved_array))

# Test fitness calculation
phenotype = zeros(3)
result = calculate_fitness!(phenotype, saved_array, opt_sys)
println("Fitness result: ", result)

# Check FFT analysis
fft_observable = saved_array[1, :]
println("Array info:")
println("  saved_array size: ", size(saved_array))
println("  fft_observable size: ", size(fft_observable))
println("  fft_observable strides: ", strides(fft_observable))
println("  time_vec length: ", length(opt_sys.t))
println("  FFT plan expecting size: ", size(opt_sys.rfft_plan))

# Try to fix the stride issue by copying
fft_observable_copy = copy(fft_observable)
println("  fft_observable_copy strides: ", strides(fft_observable_copy))

# Test what the plan was actually created for
println("\\nFFTW Plan debugging:")
plan_test_array = @view reinterpret(reshape, Float64, fill((0.0, 0.0), length(opt_sys.t)))[1, :]
println("  plan_test_array size: ", size(plan_test_array))
println("  plan_test_array strides: ", strides(plan_test_array))
println("  typeof(plan_test_array): ", typeof(plan_test_array))
println("  typeof(fft_observable_copy): ", typeof(fft_observable_copy))

# Check what type saved_array actually is
println("\\nSaved array type analysis:")
println("  typeof(saved_array): ", typeof(saved_array))
println("  typeof(saved_array[1, :]): ", typeof(saved_array[1, :]))
println("  saved_array[1, :] strides: ", strides(saved_array[1, :]))

# The issue might be that saved_array[1, :] doesn't have stride 2
# Let's try to mimic the original plan creation exactly
println("\\nTrying to match saved_array[1, :] type:")
using FFTW
if strides(saved_array[1, :]) == (2,)
    println("  Saved array has stride 2 - plan should work")
    fftData = getFrequencies(saved_array[1, :], opt_sys.rfft_plan)
else
    println("  Saved array has stride $(strides(saved_array[1, :])) - creating matching plan")
    fresh_plan = plan_rfft(saved_array[1, :])
    fftData = getFrequencies(saved_array[1, :], fresh_plan)
end
fft_peaks = find_fft_peaks(fftData)
println("FFT peaks found: ", length(fft_peaks))

if length(fft_peaks) >= 2
    println("FFT peaks at indices: ", fft_peaks)
    println("FFT peak values: ", fftData[fft_peaks])
else
    println("Insufficient FFT peaks - this explains zero fitness")
end

# Check oscillation detection
td_observable = saved_array[2, :]
is_oscillatory = check_oscillatory(td_observable)
println("Is oscillatory: ", is_oscillatory)