# Tests that the ODE solver is consistently producing the same results over time

using Test, Random
using OscillatorOptimization
using OrdinaryDiffEq: Rodas5, solve

@testset "ODE solver test" begin
    sol = solve(make_odeprob(), Rodas5(), saveat = 0.1, abstol = 1e-7, reltol = 1e-7, save_on = true)
    Amem = sol[:Amem] 
    @test argmax(Amem) == 25
end