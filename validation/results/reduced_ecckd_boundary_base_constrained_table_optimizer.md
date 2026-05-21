# Reduced ecCKD Boundary-Base Constrained Table Optimizer

Status: **constrained_table_optimizer_improved**

| Field | Value |
|---|---:|
| Base mode | boundary_weight_refit |
| Candidate scope | all_global_residual_probe |
| Residual mode | toa |
| Include Rayleigh candidates | false |
| Candidate count | 16 |
| Probe step | 0.0009765625 |
| Max log scale | 0.00390625 |
| Base objective | 7.29742558538 |
| Best exact objective | 7.2974239158 |
| Best objective reduction | 1.66957988768e-06 |
| Best TOA forcing error | 2.18922717474 W m^-2 |
| Best surface forcing error | 2.18230029605 W m^-2 |
| Accepted | true |

This artifact is separate from the canonical constrained-table optimizer.
It starts from the boundary-aware post-constrained weight-refit model
and records only boundary-base table continuations.
