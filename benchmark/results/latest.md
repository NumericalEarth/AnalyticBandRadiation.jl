# Analytic Column Benchmark

| Metric | Value |
|---|---:|
| Case | analytic_column |
| Backend | CPU |
| Precision | Float64 |
| Columns | 1 |
| Layers | 32 |
| Runtime minimum | 0.115837 ms |
| Runtime median | 0.116173 ms |
| Runtime p90 | 0.132297 ms |
| Allocations | 0 bytes |


# Staged ecCKD Cloudless Column Benchmark

| Metric | Value |
|---|---:|
| Case | staged_ecckd_cloudless_column |
| Backend | CPU |
| Precision | Float64 |
| Columns | 1 |
| Layers | 32 |
| LW g-points | 4 |
| SW g-points | 3 |
| Gas optics median | 0.001482 ms |
| Longwave solver median | 0.001945 ms |
| Shortwave solver median | 0.001117 ms |
| Heating median | 0.000082 ms |
| Total median | 0.007588 ms |
| Gas optics allocations | 0 bytes |
| Longwave solver allocations | 0 bytes |
| Shortwave solver allocations | 0 bytes |
| Heating allocations | 0 bytes |
| Total allocations | 0 bytes |


# RCEMIP-Style Column-Batch Benchmark

This is a local scaffold for the future Breeze+RRTMGP production
comparison. It runs a non-spinup, nontrivial 3D column-batch radiation
update with RCEMIP-like moist thermodynamic structure, but it does not yet
include Breeze timestep overhead, an RRTMGP baseline, GPU execution, or
Nsight profiles.

| Metric | Value |
|---|---:|
| Case | rcemip_style_column_batch |
| Status | scaffold_no_rrtmgp_baseline |
| Backend | CPU |
| Grid | 16 x 16 |
| Columns | 256 |
| Layers | 64 |
| Radiation cadence | 3600 s |
| Dynamics timestep | 300 s |
| Gas optics median | 4.523737 ms |
| Longwave solver median | 1.039119 ms |
| Shortwave solver median | 0.628199 ms |
| Heating median | 0.039361 ms |
| Radiation update median | 6.391897 ms |
| Radiation update allocations | 0 bytes |
| RRTMGP baseline | not yet measured |
| Speedup vs RRTMGP | not yet measured |
