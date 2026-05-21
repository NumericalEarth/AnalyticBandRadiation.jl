# Reduced ecCKD Retained Structural Continuation 2

Status: **constrained_table_optimizer_improved**

| Field | Value |
|---|---:|
| Base mode | boundary_table_continuation |
| Candidate scope | all_global_residual_probe |
| Residual mode | all_shortwave |
| Include Rayleigh candidates | true |
| Candidate count | 12 |
| Probe step | 0.00048828125 |
| Max log scale | 0.00390625 |
| Base objective | 7.22363619225 |
| Best ridge lambda | 0.0001 |
| Best exact objective | 7.20640276684 |
| Best objective reduction | 0.0172334254052 |
| Best TOA forcing error | 2.16192083005 W m^-2 |
| Best surface forcing error | 1.99434663304 W m^-2 |
| Accepted moves | 12 |
| Accepted | true |

This retained continuation starts from the current reduced model after
the first retained structural continuation and records another
smaller-trust-region all-shortwave residual constrained table step.
