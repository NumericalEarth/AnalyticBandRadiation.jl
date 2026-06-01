# Benchmark reports

Run from the repository root:

```sh
julia --project=. benchmark/benchmark_suite.jl
```

The command writes:

- `benchmark/results/latest.json`
- `benchmark/results/latest.md`

## Current cases

`analytic_column` is the existing single-column analytic-band timing.

`staged_ecckd_cloudless_column` times the staged runtime path for one column:
gas optics, longwave solver, shortwave solver, heating-rate conversion, and the
total pipeline. It is a CPU smoke benchmark and allocation regression check.

`rcemip_style_column_batch` is a local scaffold for the future production
Breeze+RRTMGP comparison. It runs a non-spinup 3D column-batch radiation update
with RCEMIP-like moist thermodynamic structure, a 16 x 16 horizontal grid, 64
vertical layers, and radiation cadence metadata.

## What the RCEMIP-style scaffold proves

- The staged gas-optics, solver, and heating path can run a nontrivial
  multi-column workload.
- The report includes decomposed timings and temporary-allocation counts.
- The report has explicit fields for the future RRTMGP baseline, speedup, and
  Nsight artifacts.

## What it does not prove yet

- It is not a Breeze integration benchmark.
- It does not include Breeze timestep overhead or tendency insertion.
- It does not compare against RRTMGP.
- It does not run on H100 or include Nsight Systems / Nsight Compute profiles.
- It does not establish the final 4x speedup claim.

Those missing items are required before any production performance claim.
