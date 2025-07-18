# OscTools Developer Context

## Overview
OscTools is a Julia package for optimizing oscillatory biological systems using Quality Diversity evolutionary algorithms. It specializes in finding parameter sets that produce desired oscillatory behavior in biochemical reaction networks modeled with Catalyst.jl.

## Core Architecture

### Model Integration via Catalyst.jl
- **ReactionSystem Objects**: All models are defined as Catalyst.jl `ReactionSystem` objects with tunable parameters and species
- **Automatic Conversion**: Systems are converted to `ODESystem` → `ODEProblem` for numerical integration
- **Bounded Parameters**: All tunable parameters/species have bounds metadata for optimization
- **Fixed Parameters**: Support for fixing certain parameters during optimization (e.g., `DF => 100.0`)

### Trait-Based Dispatch System (`api/traits.jl`)
**Critical for adding new models**: The trait system automatically selects fitness observables based on the reaction system name.

```julia
# Current implementations:
fitness_observables(::Val{:fullrn}) = (:Amem_old, :Amem)
fitness_observables(::Val{:trimer_rn}) = (:Amem_old, :TrimerYield)

# To add a new model:
fitness_observables(::Val{:your_model_name}) = (:fft_observable, :time_domain_observable)
```

**Requirements for new models:**
1. Both observables must be defined in the `@observables` block of your ReactionSystem
2. First observable (fft_observable) is used for FFT-based fitness scoring
3. Second observable (time_domain_observable) is used for period/amplitude calculation

### Population Generation (`api/population_generation.jl`)
- **DimArray Structure**: Uses DimensionalData.jl with dimensions `Genes` (tunable parameters/species) and `Individuals`
- **Log-Uniform Sampling**: Parameters sampled from log10 distributions within bounds
- **Constraint Handling**: Optional constraint functions during population generation
- **Fixed Parameter Exclusion**: Fixed parameters are excluded from the population matrix

### Oscillation Detection

#### FFT-Based Fitness (`api/fitness_functions/FitnessFunction.jl:1-29`)
```julia
# Fitness = standard_deviation + sum_diff
fitness = getSTD(fft_peakindexes, fftData) + getDif(fft_peakvals)
```
- Uses first observable (e.g., `Amem_old`) for FFT analysis
- Requires ≥2 peaks in frequency domain
- Combines peak variability and peak differences

#### Time-Domain Analysis (`api/fitness_functions/FitnessFunction.jl:30-51`)
- Uses second observable (e.g., `Amem`, `TrimerYield`) for period/amplitude calculation
- **Oscillatory Check**: `std(last_200_points) >= 1e-6`
- **Regularity Check**: Consecutive peaks must have consistent heights/timing within 1% variation
- **Period**: Time difference between last two peaks
- **Amplitude**: Height difference between last peak and trough

### Quality Diversity Algorithm (`evolutionary_overloads/quality-diversity.jl`)
- **Custom QD Struct**: Extends Evolutionary.jl with additional optimization system reference
- **Multi-Objective Output**: Each individual evaluated for [fitness, period, amplitude]
- **QDState**: Maintains population-wide objective values matrix (3 × population_size)
- **Parallel Evaluation**: Thread-based evaluation with progress bars
- **Constraint Integration**: Works with `WorstFitnessConstraints`

## Key Data Structures

### OptimizationReactionSystem (`api/optimization.jl:5-18`)
```julia
struct OptimizationReactionSystem{RS, OP, P, U, G, S, F, A, O}
    rx_sys::RS                # Original ReactionSystem
    oprob::OP                 # ODEProblem
    alg::A                    # ODE algorithm (e.g., Rodas5P)
    parameter_setter::P       # setp function for parameters
    species_setter::U         # setu function for species
    amem_getter::G           # getu function for observables
    amem_saver::S             # SavingCallback function
    rfft_plan::F             # Pre-planned FFT for performance
    t::Vector{Float64}        # Time vector
    fixed_params::Dict{Symbol, Float64}
    constraint_set::ConstraintSet
    obs_symbols::O            # Tuple of (fft_sym, td_sym)
end
```

### Individual Representation
- **DimVector**: Each individual is a labeled vector with `Genes` dimension
- **Metadata**: Contains selectors for parameters vs species (`u0`/`p` metadata)
- **Value Extraction**: `get_u0_p(individual)` returns `(u0=species_values, p=parameter_values)`

## File Organization

### Core API (`api/`)
- `optimization.jl`: Main optimization setup and `run_optimization`
- `evaluate_individual.jl`: ODE solving and fitness evaluation per individual
- `traits.jl`: Model-specific observable selection
- `population_generation.jl`: Initial population creation
- `constraints.jl`: Constraint definitions and application
- `results.jl`: Post-optimization data processing

### Models (`models/`)
- `full_model.jl`: Lipid oscillator models (`base_rn`, `peripheral_rn`, `fullrn`)
- `trimer_model.jl`: Trimer assembly model (`trimer_rn`)
- `accessor_functions.jl`: Catalyst.jl utilities for parameter/species extraction

### Evolutionary Overloads (`evolutionary_overloads/`)
- `quality-diversity.jl`: Custom QD algorithm implementation
- `trace.jl`: Progress tracking and visualization during optimization

## Parameter Ranges and Defaults

### Physiological Ranges (`models/full_model.jl:1-18`)
```julia
const KF_RANGE = (1e-3, 1e2)     # Forward rates (1/µM/s)
const KR_RANGE = (1e-3, 1e3)     # Reverse rates (1/s) 
const KCAT_RANGE = (1e-3, 1e3)   # Catalytic rates (1/s)
const DF_RANGE = (1.0, 1e4)      # Dimensionality factor
const L_RANGE = (1e-1, 1e2)      # Species concentrations (µM)
```

### Critical Observables
- **Amem_old**: Membrane-bound fraction (numerator excludes cytosolic complexes)
- **Amem**: Full membrane-bound fraction (includes all AP2-containing complexes)
- **TrimerYield**: Fraction of possible trimers that are assembled
- **Tmem**: Fraction of trimers that are membrane-bound

## Adding New Models

1. **Create ReactionSystem** with `@observables` block defining two observables
2. **Add trait method**: `fitness_observables(::Val{:your_model_name}) = (:obs1, :obs2)`
3. **Test integration**: Use `OptimizationReactionSystem(your_model)` constructor
4. **Verify observables**: Check that `opt_sys.obs_symbols` matches your trait

## Performance Considerations

### FFT Optimization (`api/optimization.jl:50-53`)
- Pre-planned RFFT for consistent array sizes: `plan_rfft(plan_array; flags = FFTW.PATIENT)`
- Real-valued FFT reduces computation vs complex FFT

### Threading (`evolutionary_overloads/quality-diversity.jl:192-211`)
- Use `:threadprogress` parallelization for progress bars
- Thread-safe evaluation with pre-allocated objective matrix

### Memory Management
- **SavedValues**: Efficient storage for ODE callbacks (avoids full solution storage)
- **Views**: Extensive use of array views to avoid copying
- **In-place operations**: Fitness evaluation modifies pre-allocated phenotype vectors

## Testing Structure (`test/`)
- `observable_trait_tests.jl`: Validates trait system and constructor integration
- Tests both `fullrn` and `trimer_rn` models
- Verifies observable selection and getter functionality

## Common Patterns

### Individual Evaluation Flow
1. `remake_odeprob(individual, optsys)` → Updates ODEProblem with individual's parameters
2. `solve(newprob, alg, callback=cb)` → Solves with SavingCallback 
3. `calculate_fitness!(phenotype, saved_array, optsys)` → Computes [fitness, period, amplitude]

### Result Processing
1. Extract unique individuals across all generations
2. Create DataFrame with parameter values + objectives + metadata
3. Verify all oscillatory solutions have valid fitness/period/amplitude

This architecture supports efficient optimization of biochemical oscillators while maintaining extensibility for new reaction systems through the trait-based dispatch pattern.

## Project Management Memories
- Never edit Project.toml or Manifest.toml manually, only via Pkg.jl itself
- When calling Julia from the command line, always include the `julia --project -O3` arguments