# Repository Guidelines

## Project Structure & Module Organization
- `src/OscillatorOptimization.jl`: Package entry; exports and includes.
- `src/models/`: Catalyst `ReactionSystem` definitions (`full_model.jl`, `trimer_model.jl`) and accessors.
- `src/api/`: Core APIs — traits, evaluation, fitness, population, optimization, results, constraints.
- `src/evolutionary_overloads/`: Extensions over Evolutionary.jl (quality-diversity, trace).
- `src/utils/`: Utilities (CSV/data handling, plotting helpers).
- `test/`: Unit/integration tests; `runtests.jl` orchestrates.
- `docs/`: Documenter.jl setup (`docs/make.jl`, `docs/Project.toml`).
- `examples/` and `quick_test.jl`: Runnable examples and smoke test.

## Build, Test, and Development Commands
- Setup deps: `julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'`
- Run tests: `julia --project -e 'using Pkg; Pkg.test()'`
- Smoke test: `julia --project quick_test.jl`
- Build docs: `julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'`
- Threads: `export JULIA_NUM_THREADS=4` to enable threaded evaluation.

## Coding Style & Naming Conventions
- Indentation: 4 spaces; no tabs. Keep lines reasonably short.
- Names: snake_case for functions/vars; CamelCase for types/modules; mutating functions end with `!` (e.g., `calculate_fitness!`).
- Exports: keep explicit in `src/OscillatorOptimization.jl`; add docstrings for public APIs.
- RNG: prefer deterministic seeds (`MersenneTwister`/`Xoshiro`) in examples/tests.

## Testing Guidelines
- Framework: `Test` stdlib with `@testset`/`@test`.
- Location: add tests under `test/` and include from `test/runtests.jl`.
- Determinism: set seeds; avoid long-running ODE grids. CI collects coverage and uploads to Codecov.

## Commit & Pull Request Guidelines
- Commits: imperative, concise; scope if useful.
  Examples: `add solver from DataFrame row`, `fix CI: pin SciML deps`, `update row solvers`.
- PRs must:
  - Pass CI (Julia 1.10/1.11; Ubuntu/macOS).
  - Include clear description and rationale; add benchmarks/plots if perf-related.
  - Update tests/docs when behavior or public API changes.
  - Respect `[compat]` in `Project.toml`; update CI pins if bumping SciML deps (e.g., DiffEqCallbacks, ModelingToolkit, SII).

## Security & Configuration Tips
- Data I/O: keep sample CSVs small; use `utils/datahandling.jl` helpers.
- Performance notes: use `BLAS.set_num_threads(1)` for fair comparisons; prefer `JULIA_NUM_THREADS` for parallel evaluations.

