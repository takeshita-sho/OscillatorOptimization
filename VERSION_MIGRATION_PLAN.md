# Version Migration Plan: Julia 1.11 & SciML Package Updates

## Overview

Migrate OscillatorOptimization from Julia 1.10 LTS with pinned SciML packages to Julia 1.11 with current package versions, while maintaining **identical numerical behavior** through Test-Driven Development.

## Background & Issues

### Current Version Constraints
- **Julia**: 1.10 LTS (locked due to numerical behavior changes in 1.11)
- **Pinned Packages** (locked for >1 year due to breaking API changes):
  - `DiffEqCallbacks = "3.9.1"` 
  - `ModelingToolkit = "9.41.0"`
  - `SymbolicIndexingInterface = "0.3.37"`
- **Catalyst.jl**: Indirectly locked by MTK pin, but also has breaking changes

### Observed Problems
1. **Numerical Changes**: Julia 1.10→1.11 changed ODESolution values for identical inputs
2. **API Changes**: DiffEqCallbacks and SymbolicIndexingInterface had poorly documented breaking changes
3. **Cascading Dependencies**: ModelingToolkit incompatibilities cascade through the ecosystem
4. **Model Definition Changes**: Latest Catalyst/MTK versions break fundamental model syntax

## Migration Strategy: Test-Driven Development

### Core Principle
**Lock in current behavior first, then migrate incrementally while maintaining numerical equivalence.**

### TDD Approach
1. **Establish Golden Standard**: Capture exact numerical outputs on Julia 1.10 + pinned packages
2. **Create Version Comparison Tests**: Test identical inputs produce identical outputs across versions
3. **Incremental Migration**: Update one component at a time while maintaining test suite
4. **Behavioral Validation**: Ensure all changes preserve scientific correctness

## Phase 1: Establish Numerical Baseline (Julia 1.10 LTS)

### 1.1 Capture Current Behavior
- [ ] **Expand existing tests** to cover more numerical scenarios
- [ ] **Create deterministic test cases** with fixed seeds for reproducibility
- [ ] **Save reference outputs** for comparison during migration
- [ ] **Test edge cases** that might be sensitive to version changes

### 1.2 Enhanced Test Suite
```julia
@testset "Numerical Baseline Tests" begin
    # Based on existing test/optimization_test.jl but more comprehensive
    
    @testset "ODE Solution Determinism" begin
        # Fixed individual parameters with known oscillatory behavior
        test_individual = create_test_individual()
        sol1 = solve_individual(test_individual, seed=1234)
        sol2 = solve_individual(test_individual, seed=1234)
        @test sol1.u ≈ sol2.u rtol=1e-14  # Identical for same seed
    end
    
    @testset "Fitness Function Determinism" begin
        # Fixed saved_array should produce identical fitness
        saved_array = load_reference_saved_array()
        fitness1 = calculate_fitness(saved_array)
        fitness2 = calculate_fitness(saved_array)
        @test fitness1 == fitness2
    end
    
    @testset "Population Generation Determinism" begin
        # Same seed should produce identical populations
        pop1 = generate_population(fullrn, 100, Dict(), seed=42)
        pop2 = generate_population(fullrn, 100, Dict(), seed=42)
        @test pop1 ≈ pop2 rtol=1e-14
    end
    
    @testset "End-to-End Optimization Determinism" begin
        # Small optimization run with fixed seed
        results1 = run_optimization(50, opt_sys, seed=5678)
        results2 = run_optimization(50, opt_sys, seed=5678)
        @test results1.df ≈ results2.df rtol=1e-10
    end
end
```

### 1.3 Create Reference Dataset
- [ ] **Run comprehensive test suite** on Julia 1.10 + pinned packages
- [ ] **Save numerical outputs** to version-controlled reference files
- [ ] **Include multiple model types** (fullrn, trimer_rn)
- [ ] **Cover parameter space** with diverse test cases

## Phase 2: Julia 1.11 Migration

### 2.1 Parallel Environment Setup
- [ ] **Create Julia 1.11 environment** alongside existing 1.10 setup
- [ ] **Install identical pinned package versions** on Julia 1.11
- [ ] **Run baseline tests** and identify numerical differences

### 2.2 Identify and Fix Julia Version Differences
- [ ] **Compare ODE solutions** between Julia 1.10 and 1.11 with identical packages
- [ ] **Investigate random number generation** changes between versions
- [ ] **Check floating point behavior** and compiler optimizations
- [ ] **Document specific differences** and required adjustments

### 2.3 Julia 1.11 Compatibility
```julia
@testset "Julia Version Compatibility" begin
    # Reference outputs from Julia 1.10
    reference_outputs = load_julia_1_10_references()
    
    @testset "Julia 1.11 vs 1.10 ODE Solutions" begin
        for test_case in reference_outputs.ode_cases
            julia_1_11_result = solve_ode_case(test_case)
            julia_1_10_result = test_case.reference_solution
            @test julia_1_11_result ≈ julia_1_10_result rtol=1e-12
        end
    end
end
```

## Phase 3: SciML Package Migration

### 3.1 DiffEqCallbacks Migration (3.9.1 → Current)
**Known Issues**: SavingCallback API changes

- [ ] **Identify breaking changes** in DiffEqCallbacks changelog
- [ ] **Update SavingCallback usage** in `evaluate_individual.jl`
- [ ] **Test callback behavior** matches exactly with new version
- [ ] **Validate saved values format** unchanged

```julia
@testset "DiffEqCallbacks Migration" begin
    @testset "SavingCallback Compatibility" begin
        # Test that new DiffEqCallbacks produces identical saved values
        old_saved_values = load_reference_saved_values()
        new_saved_values = generate_saved_values_new_api(same_parameters)
        @test old_saved_values ≈ new_saved_values rtol=1e-14
    end
end
```

### 3.2 SymbolicIndexingInterface Migration (0.3.37 → Current)
**Known Issues**: Parameter/species access API changes

- [ ] **Review SymbolicIndexingInterface changelog** for breaking changes
- [ ] **Update `getu`, `setp`, `setu` usage** throughout codebase
- [ ] **Test getter/setter functions** produce identical behavior
- [ ] **Validate observable access** unchanged

```julia
@testset "SymbolicIndexingInterface Migration" begin
    @testset "Parameter Setting Compatibility" begin
        # Test parameter/species setters work identically
        test_individual = create_test_individual()
        old_prob = create_ode_problem_old_api(test_individual)
        new_prob = create_ode_problem_new_api(test_individual)
        @test old_prob.p ≈ new_prob.p rtol=1e-14
        @test old_prob.u0 ≈ new_prob.u0 rtol=1e-14
    end
end
```

### 3.3 ModelingToolkit Migration (9.41.0 → Current)
**Known Issues**: Cascading incompatibilities, API changes

- [ ] **Identify specific MTK breaking changes** affecting our usage
- [ ] **Update ReactionSystem → ODESystem conversion** if needed
- [ ] **Test model compilation** produces identical symbolic systems
- [ ] **Validate parameter bounds extraction** unchanged

```julia
@testset "ModelingToolkit Migration" begin
    @testset "ReactionSystem Conversion Compatibility" begin
        # Test that new MTK produces identical ODESystem
        old_osys = convert_reaction_system_old_mtk(fullrn)
        new_osys = convert_reaction_system_new_mtk(fullrn)
        @test structural_simplify(old_osys) == structural_simplify(new_osys)
    end
end
```

## Phase 4: Catalyst.jl Migration

### 4.1 Model Definition Updates
**Known Issues**: Breaking changes in model syntax, observable definitions

- [ ] **Update `@reaction_network` syntax** to current Catalyst version
- [ ] **Fix `@observables` block syntax** if changed
- [ ] **Update parameter/species metadata** syntax
- [ ] **Test model compilation** produces equivalent behavior

### 4.2 API Compatibility
- [ ] **Update accessor functions** for parameters/species
- [ ] **Fix any changed function signatures** in Catalyst
- [ ] **Test bounds extraction** still works correctly
- [ ] **Validate tunable parameter detection** unchanged

## Phase 5: Integration Testing & Validation

### 5.1 Comprehensive End-to-End Testing
- [ ] **Run full optimization suite** on updated stack
- [ ] **Compare results with baseline** from Phase 1
- [ ] **Test performance characteristics** unchanged
- [ ] **Validate scientific correctness** of results

### 5.2 Migration Validation Framework
```julia
@testset "Migration Validation" begin
    @testset "Numerical Equivalence" begin
        # Compare old vs new stack on identical problems
        old_results = load_baseline_results()
        new_results = run_optimization_new_stack(same_parameters)
        
        @test new_results.fitness_values ≈ old_results.fitness_values rtol=1e-10
        @test new_results.periods ≈ old_results.periods rtol=1e-10
        @test new_results.amplitudes ≈ old_results.amplitudes rtol=1e-10
    end
    
    @testset "Performance Equivalence" begin
        # Ensure migration doesn't significantly impact performance
        old_benchmark = load_performance_baseline()
        new_benchmark = benchmark_new_stack()
        @test new_benchmark.time < old_benchmark.time * 1.1  # Within 10%
    end
end
```

## Implementation Approach

### Sequential Migration Strategy
1. **Never break existing functionality** - always maintain working version
2. **One package at a time** - isolate sources of breakage
3. **Extensive testing at each step** - catch issues early
4. **Document all changes** - create migration guide for future updates

### Testing Infrastructure
- **Automated CI/CD** testing on multiple Julia versions
- **Reference data version control** for numerical comparisons
- **Performance benchmarking** to catch regressions
- **Deterministic test seeds** for reproducible results

### Risk Mitigation
- **Parallel development branches** for each migration phase
- **Rollback capability** if migration introduces issues
- **Staged deployment** with thorough validation at each step
- **Community engagement** for SciML ecosystem compatibility issues

## Success Criteria

### Numerical Equivalence
- [ ] All test cases pass with `rtol=1e-12` or better
- [ ] Optimization results identical for same random seeds
- [ ] No degradation in scientific accuracy

### Performance Equivalence  
- [ ] Runtime within 10% of baseline performance
- [ ] Memory usage not significantly increased
- [ ] Compilation time acceptable

### Maintainability Improvement
- [ ] No more pinned package versions
- [ ] Compatible with latest Julia LTS
- [ ] Access to latest SciML ecosystem features
- [ ] Future-proof against dependency conflicts

## Timeline Estimate

- **Phase 1** (Baseline): 1-2 weeks
- **Phase 2** (Julia 1.11): 1 week  
- **Phase 3** (SciML packages): 2-3 weeks
- **Phase 4** (Catalyst): 1-2 weeks
- **Phase 5** (Integration): 1 week

**Total**: 6-9 weeks of focused development

## Long-term Benefits

1. **Access to latest features** in SciML ecosystem
2. **Better performance** from Julia 1.11 improvements
3. **Reduced maintenance burden** from version conflicts
4. **Community compatibility** with current package versions
5. **Future-proofing** against further ecosystem changes

This migration plan ensures that users can confidently update their Julia and package versions while maintaining the exact numerical behavior they depend on for their scientific work.