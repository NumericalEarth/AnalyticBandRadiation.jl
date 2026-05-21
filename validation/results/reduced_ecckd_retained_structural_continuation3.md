# Reduced ecCKD Retained Structural Continuation 3

Status: **constrained_table_optimizer_improved**

| Field | Value |
|---|---:|
| Base mode | boundary_table_continuation |
| Candidate scope | all_global_residual_probe |
| Residual mode | all_shortwave |
| Include Rayleigh candidates | true |
| Candidate count | 12 |
| Probe step | 0.000244140625 |
| Max log scale | 0.001953125 |
| Base objective | 7.20640276684 |
| Best ridge lambda | 0.0001 |
| Best exact objective | 7.19781262031 |
| Best objective reduction | 0.008590146534 |
| Best TOA forcing error | 2.15934378609 W m^-2 |
| Best surface forcing error | 2.01670606591 W m^-2 |
| Accepted moves | 12 |
| Accepted | true |

This retained continuation starts from the current reduced model after
the second retained structural continuation and records another
smaller-trust-region all-shortwave residual constrained table step.
