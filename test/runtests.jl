using Test
 
@testset "OscTools internal tests" begin
    include("observable_trait_tests.jl")
    include("solveODE_test.jl")
    include("optimization_test.jl")
end 