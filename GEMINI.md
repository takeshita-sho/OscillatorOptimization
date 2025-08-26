# Gemini Context for OscillatorOptimization.jl

This document provides a summary of the `OscillatorOptimization.jl` package, a Julia-based framework for optimizing oscillatory biological systems using evolutionary algorithms.

## Project Overview

`OscillatorOptimization.jl` is a specialized package developed for a computational biophysics PhD project. Its primary purpose is to discover parameter sets in biochemical reaction networks, modeled as systems of ordinary differential equations (ODEs), that produce sustained, regular oscillations.

The package is built around a custom **Quality Diversity (QD)** evolutionary algorithm. Unlike standard genetic algorithms that seek a single optimal solution, the QD approach aims to find a diverse archive of high-performing solutions, exploring various ways a system can oscillate.

## Core Architecture

The architecture is designed to be modular and extensible, particularly for incorporating new biological models.

### 1. Models as `ReactionSystem`s

-   Biological models are defined using `Catalyst.jl`'s `ReactionSystem` macro. This provides a symbolic, human-readable way to specify species, parameters, and reactions.
-   The package includes two primary built-in models:
    1.  **`fullrn` (Lipid Oscillator Model):** A comprehensive model of a lipid-based biochemical oscillator. The key output (observable) is `Amem`, representing the fraction of a specific protein (AP2) bound to the cell membrane.
    2.  **`trimer_rn` (Trimer Assembly Model):** An extension of the lipid oscillator that couples its dynamics to the assembly of a protein heterotrimer. This model introduces additional complexity and new observables like `TrimerYield`.

### 2. Trait-Based Dispatch for Model Agnosticism

-   A key feature is the use of a **trait-based dispatch system** (`src/api/traits.jl`). The `fitness_observables` trait allows the optimization pipeline to be model-agnostic.
-   By implementing this trait for a new `ReactionSystem`, a user can specify which observables should be used for fitness scoring (one for frequency-domain analysis via FFT, and one for time-domain analysis). This allows the entire framework to be applied to new models without changing the core optimization code.

### 3. The Optimization Workflow

The end-to-end workflow is as follows:

1.  **Initialization:** An `OptimizationReactionSystem` struct is created, which bundles the `ReactionSystem`, the corresponding `ODEProblem`, the ODE solver algorithm, and other metadata. The `fitness_observables` trait is called here to determine which outputs to monitor.
2.  **Population Generation:** A `DimArray` (from `DimensionalData.jl`) is used to create an initial population of individuals. Each individual is a vector of parameter values ("genes").
3.  **Evaluation Loop:**
    -   For each individual, `evaluate_individual!` (`src/api/evaluate_individual.jl`) is called.
    -   This function `remake`s the `ODEProblem` with the individual's parameters and solves it.
    -   The resulting time-series solution is passed to `calculate_fitness!` (`src/api/fitness_functions/FitnessFunction.jl`).
4.  **Fitness Calculation:**
    -   The fitness function analyzes the primary observable's time-series.
    -   It performs a Fast Fourier Transform (FFT) to identify dominant frequencies. A strong, single peak is rewarded.
    -   It also analyzes the time-domain signal to find peaks and troughs, calculating the period and amplitude of the oscillation.
    -   The final fitness score reflects the regularity and prominence of the oscillation.
5.  **Evolution:**
    -   The custom **Quality Diversity (QD) algorithm** (`src/evolutionary_overloads/quality-diversity.jl`) drives the evolution.
    -   It uses tournament selection to choose parents, and Simulated Binary Crossover (SBX) and Polynomial Mutation (PLM) to create offspring.
    -   The process is repeated over many generations to evolve a population of highly oscillatory parameter sets.
6.  **Results:** The entire run is tracked, and the final population of unique, high-fitness individuals is compiled into a `DataFrame` for analysis.

## Key Files

-   `src/OscillatorOptimization.jl`: Main module file, exports key functionalities.
-   `src/models/*.jl`: Definitions of the `ReactionSystem` models.
-   `src/api/traits.jl`: The core of the model-agnostic design.
-   `src/api/optimization.jl`: Contains the main `run_optimization` function and the `OptimizationReactionSystem` struct.
-   `src/evolutionary_overloads/quality-diversity.jl`: Implementation of the custom QD algorithm.
-   `src/api/evaluate_individual.jl`: The function that orchestrates the evaluation of a single parameter set.
-   `src/api/fitness_functions/FitnessFunction.jl`: The fitness calculation logic.
-   `Project.toml`: Lists the extensive set of dependencies, including `Catalyst`, `OrdinaryDiffEq`, `Evolutionary`, and `DimensionalData`.
