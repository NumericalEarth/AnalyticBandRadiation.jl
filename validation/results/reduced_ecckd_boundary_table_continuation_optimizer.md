# Reduced ecCKD Boundary-Table Continuation Optimizer

Status: **constrained_table_optimizer_improved**

| Field | Value |
|---|---:|
| Base mode | boundary_table_post_descent |
| Candidate scope | all_global_objective_probe |
| Residual mode | toa |
| Include Rayleigh candidates | true |
| Candidate count | 64 |
| Probe step | 0.00390625 |
| Max log scale | 0.125 |
| Base objective | 7.29585929551 |
| Best exact objective | 7.29585913012 |
| Best objective reduction | 1.65394169471e-07 |
| Best TOA forcing error | 2.18875773904 W m^-2 |
| Best surface forcing error | 2.18224163998 W m^-2 |
| Accepted | true |

This artifact starts from the retained boundary-table post-descent
state and records only additional continuation moves. It is consumed
after the boundary-base, coordinate-scan, pair-scan, and coordinate
descent artifacts.
