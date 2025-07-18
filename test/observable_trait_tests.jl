using Test
using OscillatorOptimization

models = Dict(:full => fullrn, :trimer => trimer_rn)

@testset "Trait mapping & constructor" begin
    for (lname, rn) in models
        @testset "$(lname)" begin
            fft_sym, td_sym = fitness_observables(rn)
            @test isa(fft_sym, Symbol)
            @test isa(td_sym, Symbol)

            # Quick constructor to ensure pipeline builds
            opt = OptimizationReactionSystem(rn; fixed_params = Dict(:DF => 100.0), tspan=(0.0, 1.0), dt=0.5)
            @test opt.obs_symbols == (fft_sym, td_sym)
            vals = opt.amem_getter(opt.oprob.u0)
            @test length(vals) == 2
        end
    end
end 