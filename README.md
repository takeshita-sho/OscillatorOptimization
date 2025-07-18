# **Geometrically Tunable Oscillator: Optimization of Oscillatory Systems Using Evolutionary Algorithms**

## **Overview**

This module implements a custom Quality Diversity (QD) evolutionary algorithm to optimize oscillatory biological systems modeled with differential equations. The project supports two main models: a lipid oscillator system and a coupled trimer assembly system. The goal is to find parameter sets that produce desired oscillatory behavior in systems of ordinary differential equations (ODEs). The codebase integrates with the [Evolutionary.jl](https://github.com/wildart/Evolutionary.jl) framework, uses trait-based dispatch for model-agnostic optimization, and leverages Julia's capabilities for high-performance numerical computing.

## **Table of Contents**

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Getting Started](#getting-started)
- [Code Structure](#code-structure)
  - [1. Reaction System Model Definition](#1-reaction-system-model-definition)
  - [2. Accessor Functions for Reaction System](#2-accessor-functions-for-reaction-system)
  - [3. Trait-Based Dispatch System](#3-trait-based-dispatch-system)
  - [4. Initial Population Generation](#4-initial-population-generation)
  - [5. Optimization Setup](#5-optimization-setup)
  - [6. Custom Quality Diversity Algorithm](#6-custom-quality-diversity-algorithm)
  - [7. Individual Evaluation](#7-individual-evaluation)
  - [8. Fitness Function](#8-fitness-function)
  - [9. Fitness Function Helpers](#9-fitness-function-helpers)
  - [10. Tracing and Recording](#10-tracing-and-recording)
  - [11. Results Processing](#11-results-processing)
- [Workflow](#workflow)
- [Usage](#usage)
- [Examples](#examples)

## **Features**

- Custom implementation of a Quality Diversity evolutionary algorithm.
- Integration with Evolutionary.jl for optimization routines.
- **Trait-based dispatch system** for model-agnostic optimization pipeline.
- **Multiple model support**: Lipid oscillator and trimer assembly models.
- Evaluation of ODE systems with customizable parameters and initial conditions.
- **Model-specific fitness observables** automatically selected via trait system.
- Fitness functions designed to assess oscillatory behavior.
- Parallel evaluation of individuals with multithreading support.
- Comprehensive tracing and recording of optimization progress.
- Results aggregation and data processing for analysis and visualization.

## **Code Structure**

The codebase is organized into several modules, each handling different aspects of the optimization process. Below is a detailed explanation of each component.

### **1. Reaction System Model Definition**

**Files:** `full_model.jl`, `trimer_model.jl`

**Purpose:**

- Defines biochemical oscillator models as `ReactionSystem` objects using Catalyst.jl, which allows for easy construction of chemical reaction networks.
- Specifies the components of the reaction systems, including rate constants, dimensionality factors (DF), species, and initial conditions.
- **Lipid Oscillator Models** (`full_model.jl`):
  - `base_rn`: 12-variable restricted model with essential complexes
  - `peripheral_rn`: Additional peripheral reactions
  - `fullrn`: Complete reaction network extending `base_rn` with `peripheral_rn`
- **Trimer Assembly Model** (`trimer_model.jl`):
  - `trimer_rn`: Extended model with heterotrimer assembly on membrane
  - Includes monomer species (B, C, D) and assembly intermediates
  - Models pathway from monomers → dimers → trimers with membrane localization

**Key Components:**

- **Lipid Models**: `base_rn`, `peripheral_rn`, `fullrn` with observables `Amem_old`, `Amem`
- **Trimer Model**: `trimer_rn` with observables `Amem_old`, `TrimerYield`, `Tmem`
- **make_odeprob/make_trimer_odeprob**: Functions that convert reaction networks into `ODEProblem` objects

**Workflow Integration:**

- Models are converted to `ODESystem` objects which form the basis for optimization. Parameters and species initial conditions are tunable, enabling evolutionary search for oscillatory behavior. The trait system automatically selects appropriate observables for each model type.

### **2. Accessor Functions for Reaction System**

**File:** `accessor_functions.jl`

**Purpose:**

- Provides utility functions to access and manipulate the symbolic properties of `ReactionSystem` objects.
- Enables easy extraction of tunable parameters and species, making the optimization process more manageable.
- Works with both lipid oscillator and trimer assembly models.

**Key Components:**

- **get_parameter_symbols:** Extracts the parameters of the reaction system, with an option to filter only tunable ones.
- **get_species_symbols:** Extracts the species of the reaction system, with an option to filter only tunable ones.
- **get_tunable_bounds_dictionary:** Returns bounds for tunable parameters in a dictionary format.
- **get_unfixed_tunable_symbols:** Returns tunable symbols excluding fixed parameters.
- **get_default_values:** Returns default values for tunable parameters.

**Workflow Integration:**

- These functions are used throughout the optimization process to access and manipulate the properties of reaction systems dynamically, especially during population initialization and evaluation. They work seamlessly with both model types through the shared parameter structure.

### **3. Trait-Based Dispatch System**

**File:** `traits.jl`

**Purpose:**

- Implements a trait-based dispatch system for model-agnostic optimization.
- Automatically selects appropriate fitness observables based on the reaction system's name.
- Enables seamless integration of new models without modifying the optimization pipeline.

**Key Components:**

- **fitness_observables trait:** Returns a tuple `(fft_observable, time_domain_observable)` for each model type.
- **Value-based dispatch:** Uses `Val{:model_name}` to dispatch on reaction system names.
- **Built-in model support:**
  - `fullrn`: Uses `(:Amem_old, :Amem)` observables
  - `trimer_rn`: Uses `(:Amem_old, :TrimerYield)` observables

**Workflow Integration:**

- The trait system is automatically invoked during `OptimizationReactionSystem` construction, ensuring the correct observables are used for fitness evaluation without requiring model-specific code changes in the optimization pipeline.

### **4. Initial Population Generation**

**File:** `population_generation.jl`

**Purpose:**

- Generates the initial population of potential solutions for the optimization.
- Uses a `DimArray` structure to represent individuals, allowing for efficient, labeled indexing.

**Key Components:**

- **generate_population:** Generates an initial population, excluding fixed inputs, and ensures each individual satisfies constraints if provided.
- **make_tunable_distributions:** Creates a distribution for each tunable parameter, allowing for sampling during population generation.
- **DimArray Representation:** Uses `DimensionalData.jl` to create labeled data structures for easy access and manipulation of individuals.

**Workflow Integration:**

- This population serves as the starting point for the evolutionary optimization. It undergoes selection, mutation, and crossover to evolve towards the desired oscillatory behavior.

### **4. Optimization Setup**

**File:** `optimization.jl`

**Purpose:**

- Initializes the optimization problem for any supported model type.
- Defines the `ObjectiveFunction` that wraps the evaluation of individuals.
- Sets up the `OptimizationReactionSystem`, which includes the ODE problem, parameter and species setters, and observable getters.
- Configures the FFT plans for efficient frequency analysis.
- **Uses trait-based dispatch** to automatically select appropriate fitness observables for each model.

**Key Components:**

- **OptimizationReactionSystem:** A struct that holds the ODE problem and related functions, with model-specific observables determined by traits.
- **ObjectiveFunction:** Wraps the evaluation function to be used by the evolutionary algorithm.
- **FFT Plans:** Pre-planned FFT computations for performance optimization.
- **Trait Integration:** Automatically determines fitness observables via `fitness_observables(reaction_system)` trait.

**Workflow Integration:**

- The constructor automatically detects the model type and selects appropriate observables (e.g., `Amem_old`/`Amem` for lipid models, `Amem_old`/`TrimerYield` for trimer models) without requiring model-specific code changes.

### **6. Custom Quality Diversity Algorithm**

**File:** `quality-diversity.jl`

**Purpose:**

- Implements a custom QD algorithm by extending and overloading methods from Evolutionary.jl.
- Defines custom structures and functions to control the evolutionary process.

**Key Components:**

- **QD Struct:** Encapsulates algorithm parameters like population size, crossover rate, mutation rate, selection, crossover, mutation functions, and metrics.
- **QDState Struct:** Maintains the state of the optimization, including the best individual and objective values.
- **Overloaded Methods:** Custom implementations of Evolutionary.jl methods such as `initial_state`, `update_state!`, `recombine!`, `mutate!`, and `evaluate!`.
- **Objective Function Integration:** Interfaces with `EvolutionaryObjective` to evaluate individuals using the custom objective function.

### **7. Individual Evaluation**

**File:** `evaluate_individual.jl`

**Purpose:**

- Defines how each individual (a set of parameters and initial conditions) is evaluated.
- Solves the ODEProblem with the individual's parameters.
- Calculates the fitness based on the solution.

**Key Components:**

- **evaluate_individual!:** Main function that evaluates an individual and updates the phenotype vector with fitness, period, and amplitude.
- **solve_odes:** Solves the ODEProblem for a given individual.
- **remake_odeprob:** Remakes the ODEProblem with new parameters and species initial conditions.

### **8. Fitness Function**

**File:** `FitnessFunction.jl`

**Purpose:**

- Calculates the fitness of an individual based on the ODE solution.
- Analyzes the time-series data to determine oscillatory behavior.

**Key Components:**

- **calculate_fitness!:** Core function that computes fitness, period, and amplitude from the time-series data.
- **find_fft_peaks:** Identifies significant peaks in the frequency domain.
- **find_amem_peaks:** Identifies peaks and troughs in the time domain.

### **9. Fitness Function Helpers**

**File:** `fitness_function_helpers.jl`

**Purpose:**

- Provides helper functions used in the fitness calculation.
- Handles signal processing tasks like FFT computation and peak analysis.

**Key Components:**

- **getDif:** Calculates the sum of absolute differences between consecutive peak values.
- **getSTD:** Computes the summed standard deviation around peaks in the FFT data.
- **getFrequencies:** Performs FFT on the time-series data.
- **compute_period:** Calculates the period between the last two peaks.
- **compute_amplitude:** Calculates the amplitude between the last peak and last trough.

### **10. Tracing and Recording**

**File:** `trace.jl`

**Purpose:**

- Records the state of the optimization at each iteration.
- Saves information about the population, fitness values, and solutions.

**Key Components:**

- **trace!:** Custom function that records various aspects of the optimization state, including populations and fittest individuals.
- **insert_fixed_inputs!:** Inserts fixed parameters back into an individual's chromosome for reporting purposes.

### **11. Results Processing**

**File:** `results.jl`

**Purpose:**

- Processes and aggregates the results after the optimization run.
- Compiles the data into a `DataFrame` for analysis.

**Key Components:**

- **Results Struct:** Holds the hit rate, results `DataFrame`, and combined population.
- **Results Function:** Processes the `EvolutionaryOptimizationResults` to construct a `Results` object, extracting unique individuals and compiling metrics.

## **Workflow**

1. **Initialization:**

   - Define the ODE model and problem using `OptimizationReactionSystem`.
   - **Trait system automatically selects** appropriate fitness observables for the model type.
   - Set up the initial population with parameters and species initial conditions.

2. **Evaluation:**

   - For each individual, `evaluate_individual!` solves the ODEs and calculates the fitness using `calculate_fitness!`.
   - Fitness calculation involves signal processing of the ODE solution to assess oscillatory behavior.
   - **Model-specific observables** are used automatically (e.g., `Amem_old`/`Amem` for lipid models, `Amem_old`/`TrimerYield` for trimer models).

3. **Evolutionary Loop:**

   - The custom QD algorithm evolves the population through selection, crossover, and mutation.
   - Overloaded methods from Evolutionary.jl control the evolutionary process.

4. **Tracing:**

   - At each iteration, `trace!` records the optimization state, including populations and fittest individuals.

5. **Results Processing:**

   - After optimization, `Results` processes the collected data.
   - Unique individuals are extracted, and results are compiled into a `DataFrame`.

6. **Analysis:**

   - The results can be analyzed and visualized to interpret the optimization outcomes.

## **Usage**

1. **Set Up the ODE Model:**

   - Define your ODE system (choose from `fullrn`, `trimer_rn`, or define your own).
   - Create an `OptimizationReactionSystem` instance - **traits automatically handle observable selection**.
   - **For new models:** Implement the `fitness_observables(::Val{:your_model_name})` trait method.

2. **Configure the Optimization:**

   - Initialize the `QD` struct with desired algorithm parameters.
   - Set up the initial population.

3. **Run the Optimization:**

   ```julia
   #- Activate the project
   using DrWatson
   @quickactivate "GeometricallyTunableOscillator"
   include(srcdir("OscTools", "OscTools.jl"))
   using .OscTools

   begin
       # Set BLAS threads to 1 to avoid multithreading issues
       using LinearAlgebra
       BLAS.set_num_threads(1)

       using OrdinaryDiffEq
       using CSV
   end

   # Define fixed parameters
   fixed_params = Dict(:DF => 100.0)

   # Define ODE solver algorithm
   alg = Rodas5P(chunk_size = 16)
   # Define optimization system for lipid oscillator model
   fixed_opt_sys = OptimizationReactionSystem(fullrn, alg; fixed_params)

   # Run optimization
   fixed_results = run_optimization(1000, fixed_opt_sys; mutation_δ = 1.0, pm = 0.75, η = 1, n_points = Inf, num_tournament_groups = 20,  mutationRate = 0.95, crossoverRate = 0.75, sbx_pm = 0.3, sbx_η = 2, parallelization = :threadprogress, callback = nothing, show_trace = true)

   # Access results dataframe and save to disk
   results_df = fixed_results.results_df
   CSV.write(datadir("optimization_results", "fixed_DF100_results.csv"), results_df)

   # Alternative: Use trimer assembly model
   # trimer_opt_sys = OptimizationReactionSystem(trimer_rn, alg; fixed_params)
   # trimer_results = run_optimization(1000, trimer_opt_sys; <same parameters>)