# OscillatorOptimization Package Generalization Plan

## Overview

Transform OscillatorOptimization from a specialized package for specific biochemical models into a generic Quality Diversity optimizer for any Catalyst.jl ReactionSystem with oscillatory behavior.

## Development Approach: Test-Driven Development (TDD)

We will use Test-Driven Development to ensure:
1. **Backward compatibility** - Current hardcoded models continue to work exactly as before
2. **Regression prevention** - No behavioral changes during refactoring
3. **User-focused design** - API designed around actual use cases
4. **Confidence in changes** - Every modification validated by tests

### TDD Process
1. **Define end-user use cases** and desired API
2. **Write tests** for both existing behavior and new generic behavior
3. **Implement changes** to make tests pass
4. **Refactor** while keeping all tests green
5. **Repeat** for each component

## Current Hardcoded Assumptions to Generalize

### 1. Hardcoded Observable Names
**Current Problem:**
```julia
# src/api/traits.jl:28-29
fitness_observables(::Val{:fullrn}) = (:Amem_old, :Amem)
fitness_observables(::Val{:trimer_rn}) = (:Amem_old, :TrimerYield)
```
Only supports specific models with predefined observable names.

**Target API:**
```julia
OptimizationReactionSystem(rx_sys; 
    fft_observable=:MyFFTObs, 
    time_domain_observable=:MyTimeObs)
```

### 2. Hardcoded Model-Specific Logic
**Current Problem:**
```julia
# src/models/accessor_functions.jl:115-122
function get_amem_symbols(rs::AbstractSystem)
    numerator_symbols = Symbol.(filter(x -> occursin("LpA", x), symbol_strings))
    denominator_symbols = Symbol.(filter(x -> occursin("A", x), symbol_strings))
```
Assumes all models have "LpA" and "A" patterns in variable names.

**Target:** Remove model-specific symbol extraction, rely on user-specified observables.

### 3. Fixed Parameter Ranges
**Current Problem:**
```julia
# src/models/full_model.jl:1-18
const KF_RANGE = (1e-3, 1e2) #forward rate constants
const BOUNDS = dictionary([:kfᴸᴬ => KF_RANGE, :krᴸᴬ => KR_RANGE, ...])
```
Only supports specific biochemical parameter types and ranges.

**Target:** Extract bounds from Catalyst metadata:
```julia
@parameters k1::Float64 = 1.0, [bounds=(0.1, 10.0)]
@species X(t) = 5.0, [bounds=(1.0, 100.0)]
```

### 4. Fixed Time Domain Settings
**Current Problem:**
```julia
# api/optimization.jl:36, 45
tspan = (0.0, 1968.2), dt = 0.1
```
Hardcoded time settings may not suit all oscillatory systems.

**Target API:**
```julia
OptimizationReactionSystem(rx_sys; 
    tspan=(0.0, 100.0), 
    dt=0.05,
    abstol=1e-8, 
    reltol=1e-6)
```

### 5. Hardcoded Fitness Function Logic
**Current Problem:**
```julia
# src/api/fitness_functions/FitnessFunction.jl:92-94
function check_oscillatory(amem)
    last_200_points = @view amem[end-2000:end]  # Assumes dt=0.1!
    return std(last_200_points) >= 1e-6
```
Hardcoded thresholds and time windows.

**Target:** Configurable oscillation detection:
```julia
struct OscillationConfig
    oscillatory_threshold::Float64     # Default: 1e-6
    regularity_threshold::Float64      # Default: 0.01
    min_analysis_points::Int           # Default: 2000
    min_peaks::Int                     # Default: 2
end
```

## Proposed Generic API Design

### Core Constructor
```julia
OptimizationReactionSystem(
    rx_sys::ReactionSystem;
    fft_observable::Symbol,
    time_domain_observable::Symbol,
    time_config::TimeConfig = TimeConfig(),
    oscillation_config::OscillationConfig = OscillationConfig(),
    fitness_function::AbstractFitnessFunction = FFTAmplitudeFitness(),
    fixed_params::Dict{Symbol, Float64} = Dict{Symbol, Float64}(),
    constraint_set::ConstraintSet = ConstraintSet([null_constraint()]),
    alg = Rodas5P(autodiff = AutoForwardDiff(chunksize=length(species(rx_sys)))),
    remove_conserved::Bool = false
)
```

### Configuration Structs
```julia
struct TimeConfig
    tspan::Tuple{Float64, Float64}
    dt::Float64
    abstol::Float64
    reltol::Float64
end

struct OscillationConfig
    oscillatory_threshold::Float64
    regularity_threshold::Float64
    min_analysis_points::Int
    min_peaks::Int
end

abstract type AbstractFitnessFunction end

struct FFTAmplitudeFitness <: AbstractFitnessFunction
    config::OscillationConfig
end
```

### User-Friendly Convenience Functions
```julia
# For backward compatibility and common use cases
function create_lipid_oscillator_optimizer(model_variant::Symbol = :fullrn)
    rx_sys = get_builtin_model(model_variant)
    fft_obs, td_obs = get_default_observables(model_variant)
    return OptimizationReactionSystem(rx_sys; 
        fft_observable=fft_obs, 
        time_domain_observable=td_obs)
end
```

## Implementation Phases

### Phase 1: Test Infrastructure
- [ ] **Set up comprehensive test suite** for existing behavior
- [ ] **Test current models** (`fullrn`, `trimer_rn`) produce identical results
- [ ] **Benchmark current performance** as regression baseline
- [ ] **Test edge cases** in current implementation

### Phase 2: Observable Parameterization
- [ ] **Add observable parameters** to `OptimizationReactionSystem` constructor
- [ ] **Remove hardcoded traits system** while maintaining backward compatibility
- [ ] **Validate observables exist** in provided ReactionSystem
- [ ] **Test equivalence** with original trait-based approach

### Phase 3: Time Configuration
- [ ] **Parameterize time settings** (tspan, dt, tolerances)
- [ ] **Update oscillation detection** to use configurable time windows
- [ ] **Add TimeConfig struct** with sensible defaults
- [ ] **Test time-dependent calculations** remain correct

### Phase 4: Bounds Generalization
- [ ] **Extract bounds from Catalyst metadata** instead of hardcoded dictionaries
- [ ] **Support arbitrary parameter/species names**
- [ ] **Validate bounds exist** for all tunable parameters
- [ ] **Test population generation** with new bounds system

### Phase 5: Fitness Function Modularity
- [ ] **Create AbstractFitnessFunction interface**
- [ ] **Refactor existing fitness logic** into FFTAmplitudeFitness
- [ ] **Make oscillation detection configurable**
- [ ] **Test fitness calculations** produce identical results

### Phase 6: User Experience & Documentation
- [ ] **Add convenience constructors** for common use cases
- [ ] **Create comprehensive examples** for different model types
- [ ] **Write tutorial documentation** for new users
- [ ] **Add error messages** to guide users when things go wrong

### Phase 7: Advanced Features
- [ ] **Support custom constraint functions**
- [ ] **Add alternative fitness functions** (frequency domain, time domain, etc.)
- [ ] **Enable multi-objective optimization** beyond period/amplitude
- [ ] **Performance optimization** for large models

## Test Strategy

### Regression Tests
- Current models (`fullrn`, `trimer_rn`) must produce **identical numerical results**
- All existing examples must continue to work without modification
- Performance must not degrade significantly

### Integration Tests
- Generic API works with user-defined ReactionSystems
- Bounds extraction works with various Catalyst models
- Observable validation catches invalid inputs early

### Unit Tests
- Each configuration struct validates inputs correctly
- Fitness functions are modular and composable
- Error messages are helpful and actionable

### Property-Based Tests
- Population generation respects bounds for any valid ReactionSystem
- Oscillation detection is consistent across different time configurations
- Fitness calculations are numerically stable

## Success Criteria

1. **Zero breaking changes** for existing users
2. **Simple API** for new users with custom models
3. **Clear documentation** with examples for different domains
4. **Extensible architecture** for future fitness functions
5. **Performance parity** with current implementation
6. **Comprehensive test coverage** (>90%)

## Example Target Usage

### For Existing Users (No Changes Required)
```julia
# Still works exactly as before
opt_sys = OptimizationReactionSystem(fullrn)
results = run_optimization(1000, opt_sys)
```

### For New Users with Custom Models
```julia
# User defines their own oscillatory model
@reaction_network my_oscillator begin
    @parameters k1 = 1.0, [bounds=(0.1, 10.0)]
                k2 = 2.0, [bounds=(1.0, 20.0)]
    @species A(t) = 5.0, [bounds=(1.0, 100.0)]
             B(t) = 2.0, [bounds=(0.1, 50.0)]
    
    @observables begin
        A_ratio ~ A / (A + B)
        total ~ A + B
    end
    
    k1, A + B --> 2A
    k2, 2A --> A + B
end

# Configure optimizer for their specific observables
opt_sys = OptimizationReactionSystem(my_oscillator;
    fft_observable=:A_ratio,
    time_domain_observable=:total,
    tspan=(0.0, 100.0),
    dt=0.1
)

results = run_optimization(500, opt_sys)
```

This plan ensures we maintain the robust, tested behavior users depend on while opening the package to entirely new domains and use cases.