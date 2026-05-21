# Reduced ecCKD Constrained Table Optimizer

Status: **constrained_table_optimizer_rejected**

| Field | Value |
|---|---:|
| Target case | ecckd_clear_sky_tropical_column |
| Target metric | toa_forcing_max_abs |
| Candidate scope | all_global_residual_probe |
| Residual mode | boundary |
| Include Rayleigh candidates | true |
| Candidate count | 96 |
| Probe step | 0.001953125 |
| Max log scale | 0.015625 |
| Base objective | 7.13985029534 |
| Best ridge lambda | 1000000 |
| Best exact objective | 7.22878936628 |
| Best objective reduction | -0.0889390709339 |
| Best TOA forcing error | 2.16652518014 W m^-2 |
| Best surface forcing error | 2.16863680988 W m^-2 |
| Accepted | false |

This is the first constrained multi-parameter table optimizer preflight.
It linearizes hard-gate-scaled shortwave flux, heating-rate, and boundary
residuals with respect to multiple active table entries, solves bounded
ridge updates, then evaluates the nonlinear full hard objective before
accepting any move.
