# Benchmarking

Independent H100 / CPU benchmark for `Lightflux.jl` against
RRTMGP.jl, driven through Breeze's `update_radiation!` so the comparison runs
under the same call surface production uses. The project path-deps onto a
developing Breeze checkout so we can evolve the gas-optics + transport kernels
in `BreezeRadiativeHeatingExt` and benchmark them in one step.

Memory rule for this benchmark: **no array of shape `Nx * Ny * Nz * Nk`**.
G-points stream inside the column transport kernel
(`BreezeRadiativeHeatingExt._tabulated_ecckd_streaming_radiation!`); optical
depth is computed inline per g-point and discarded after the up/down sweeps,
so persistent storage scales as `Nx * Ny * Nz` only.

## Layout

- `Project.toml` — environment. `[sources]` path-deps onto this repo's
  `Lightflux.jl` and the developing Breeze checkout at
  `/shared/home/greg/Projects/BreezeRadiativeHeatingDev/Breeze.jl`.
- `rcemip_h100_benchmark.jl` — the benchmark script. ENV-configurable for
  k-model, grid, samples, output dir.
- `h100_kmodel_coverage.sbatch` — Slurm job that runs the benchmark for five
  validated ecCKD k-model combinations at the same 512×512×128 / 262,144
  column workload.
- `results/` — per-run output, one subdirectory per `rcemip_h100_<label>`.

## Quick start

```sh
julia --project=benchmarking -e 'using Pkg; Pkg.instantiate()'

# small CPU smoke
RADIATIVE_HEATING_BACKEND=CPU \
RADIATIVE_HEATING_NX=4 RADIATIVE_HEATING_NY=4 RADIATIVE_HEATING_NZ=16 \
RADIATIVE_HEATING_SAMPLES=2 \
RADIATIVE_HEATING_GAS_MODEL_SOURCE=validated_ecCKD \
RADIATIVE_HEATING_ECCKD_LW_PATH=$(julia --project=. -e 'using Lightflux; print(official_ecckd_definition_path(:longwave_32))') \
RADIATIVE_HEATING_ECCKD_SW_PATH=$(julia --project=. -e 'using Lightflux; print(official_ecckd_definition_path(:shortwave_32))') \
RADIATIVE_HEATING_OUTPUT_DIR=benchmarking/results/cpu_smoke \
julia --project=benchmarking benchmarking/rcemip_h100_benchmark.jl

# H100 production sweep
sbatch benchmarking/h100_kmodel_coverage.sbatch
```

## Environment variables

| Variable                                | Default                  | Meaning                                            |
|----------------------------------------|--------------------------|----------------------------------------------------|
| `RADIATIVE_HEATING_BACKEND`            | `CPU`                    | `CPU` or `H100`                                    |
| `RADIATIVE_HEATING_GAS_MODEL_SOURCE`   | `validated_ecCKD`        | `validated_ecCKD` or `synthetic_fixed_coefficients`|
| `RADIATIVE_HEATING_ECCKD_LW_PATH`      | (required for validated) | absolute path to LW `ecckd-*-definition.nc`        |
| `RADIATIVE_HEATING_ECCKD_SW_PATH`      | (required for validated) | absolute path to SW `ecckd-*-definition.nc`        |
| `RADIATIVE_HEATING_NG_LW` / `NG_SW`    | `32` / `16`              | g-points for synthetic-coeff path only             |
| `RADIATIVE_HEATING_NX` / `NY` / `NZ`   | `32` / `32` / `64`       | grid dimensions                                    |
| `RADIATIVE_HEATING_SAMPLES`            | `5`                      | post-warmup timing samples                         |
| `RADIATIVE_HEATING_RUN_RRTMGP`         | `true`                   | also run the RRTMGP baseline                       |
| `RADIATIVE_HEATING_OUTPUT_DIR`         | `benchmarking/results/…` | output directory for the JSON/MD artifacts         |
| `RADIATIVE_HEATING_LABEL`              | `""`                     | optional suffix on the auto-generated output dir   |
| `RADIATIVE_HEATING_CO2_VMR` (and friends) | published defaults     | trace-gas mole fractions used at every column      |

## Output

Each run writes `radiative_heating_rcemip_latest.{json,md}` to its output
directory. `figures/make_pr_figures.jl` reads those JSONs to populate the H100
speedup figure in the PR.
